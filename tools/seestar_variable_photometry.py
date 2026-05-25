#!/usr/bin/env python3
"""
Seestar S30 Pro — Variable Star Photometry Pipeline
====================================================
Génère des courbes de lumière différentielles à partir des FITS bruts
du Seestar en utilisant Siril comme moteur de traitement.

Workflow complet :
    1. Lien symbolique des frames FITS du Seestar → séquence Siril
    2. Calibration optionnelle (darks / flats)
    3. Résolution astrométrique de chaque frame (seqplatesolve)
    4. Alignement de la séquence (seqapplyreg)
    5. Photométrie d'ouverture différentielle (light_curve)
    6. Export CSV compatible AAVSO

Prérequis :
    pip install astropy astroquery
    Siril 1.4+ installé avec siril-cli dans le PATH

Usage rapide :
    python seestar_variable_photometry.py --lights /data/seestar/lights/
    python seestar_variable_photometry.py --lights /data/seestar/lights/ \\
        --target "RR Lyr" --refs "19:25:27.9,+42:47:04" "19:24:50.2,+42:50:12"
    python seestar_variable_photometry.py --lights /data/seestar/lights/ \\
        --config my_star.json          # Reprendre une configuration sauvegardée
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ── Dépendances optionnelles ──────────────────────────────────────────────────
try:
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    HAS_ASTROPY = True
except ImportError:
    HAS_ASTROPY = False

try:
    from astroquery.simbad import Simbad
    HAS_SIMBAD = True
except ImportError:
    HAS_SIMBAD = False


# ── Profil matériel Seestar S30 Pro ───────────────────────────────────────────
SEESTAR_S30_PRO = {
    "focal_mm": 150,          # Longueur focale en mm
    "pixel_um": 2.9,          # Taille pixel en µm
    "gain_e_adu": 80,         # Gain électrons/ADU (approximatif)
    "oscsensor": "ZWO Seestar S30",
    "oscfilter_bb": "UV/IR Block",
    "oscfilter_nb": "ZWO Seestar LP",
    "bayer": "GRBG",
    # FOV approximatif : atan2(sensor_size/focal) ≈ 1.5° × 1.1°
    "fov_deg": 1.5,
}

# ── Paramètres photométrie par défaut ─────────────────────────────────────────
PHOT_DEFAULTS = {
    # Canal couleur pour OSC : 1 = Vert (meilleur SNR sur capteur Bayer GRBG)
    # 0=R  1=G  2=B  -1=luminance combinée
    "channel": 1,
    "aperture_px": 10,   # Rayon d'ouverture en pixels
    "inner_px": 15,      # Rayon interne de l'anneau de ciel
    "outer_px": 25,      # Rayon externe de l'anneau de ciel
}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("seestar_phot")


# ─────────────────────────────────────────────────────────────────────────────
# Siril runner : génère et exécute des scripts .ssf via siril-cli
# ─────────────────────────────────────────────────────────────────────────────

def find_siril_cli() -> str:
    """Retourne le chemin de siril-cli ou lève FileNotFoundError."""
    candidates = [
        "siril-cli",
        "/usr/bin/siril-cli",
        "/usr/local/bin/siril-cli",
        "/opt/homebrew/bin/siril-cli",
        "/Applications/Siril.app/Contents/MacOS/siril-cli",
    ]
    for c in candidates:
        if shutil.which(c):
            return c
    raise FileNotFoundError(
        "siril-cli introuvable.\n"
        "  → Installez Siril 1.4+ depuis https://siril.org\n"
        "  → Vérifiez que siril-cli est dans votre PATH"
    )


def run_siril_script(working_dir: Path, commands: List[str], script_name: str = "_pipeline.ssf") -> bool:
    """
    Écrit un script .ssf et l'exécute via siril-cli.
    Retourne True si succès, False si erreur.
    """
    script_path = working_dir / script_name
    content = "requires 1.4\n" + "\n".join(commands) + "\n"
    script_path.write_text(content, encoding="utf-8")

    log.debug("Script Siril :\n%s", content)

    try:
        cli = find_siril_cli()
    except FileNotFoundError as e:
        log.error("%s", e)
        return False

    cmd = [cli, "-d", str(working_dir), "-s", str(script_path)]
    log.info("Exécution Siril : %s", " ".join(cmd))

    proc = subprocess.run(cmd, capture_output=True, text=True)

    for line in proc.stdout.splitlines():
        log.info("[siril] %s", line)
    for line in proc.stderr.splitlines():
        level = log.error if "error" in line.lower() else log.debug
        level("[siril] %s", line)

    if proc.returncode != 0:
        log.error("siril-cli a retourné le code %d", proc.returncode)
        return False
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Utilitaires coordonnées
# ─────────────────────────────────────────────────────────────────────────────

def resolve_star_name(name: str) -> Optional[Tuple[float, float]]:
    """Résout un nom d'étoile → (RA_deg, Dec_deg) via Simbad."""
    if not HAS_SIMBAD:
        log.warning("astroquery non installé — résolution Simbad indisponible. pip install astroquery")
        return None
    try:
        log.info("Recherche '%s' dans Simbad...", name)
        result = Simbad.query_object(name)
        if result is None:
            return None
        coord = SkyCoord(result["RA"][0], result["DEC"][0], unit=(u.hourangle, u.deg))
        return coord.ra.deg, coord.dec.deg
    except Exception as exc:
        log.warning("Simbad: %s", exc)
        return None


def parse_radec(value: str) -> Tuple[float, float]:
    """
    Convertit une chaîne 'RA,Dec' en (RA_deg, Dec_deg).
    Formats acceptés :
      - Degrés décimaux  : "83.8221,-5.3911"
      - Sexagésimal      : "05:35:17.3,-05:23:28"
    """
    parts = value.split(",")
    if len(parts) != 2:
        raise ValueError(f"Format attendu 'RA,Dec', reçu : {value!r}")
    ra_s, dec_s = parts[0].strip(), parts[1].strip()
    if HAS_ASTROPY:
        if ":" in ra_s:
            coord = SkyCoord(ra_s, dec_s, unit=(u.hourangle, u.deg))
        else:
            coord = SkyCoord(float(ra_s), float(dec_s), unit=u.deg)
        return coord.ra.deg, coord.dec.deg
    return float(ra_s), float(dec_s)


def wcs_arg(ra: float, dec: float) -> str:
    """Formate RA/Dec pour l'argument Siril -wcs= ou -refwcs="""
    return f"{ra:.6f},{dec:+.6f}"


# ─────────────────────────────────────────────────────────────────────────────
# Pipeline principal
# ─────────────────────────────────────────────────────────────────────────────

class PhotometryPipeline:
    """
    Pipeline de photométrie pour étoiles variables avec le Seestar S30 Pro.

    Structure des répertoires créés :
        working_dir/
        ├── process/          Frames liées et traitées
        ├── masters/          Darks/flats masterisés
        └── results/          light_curve.dat, PNG, CSV AAVSO
    """

    def __init__(
        self,
        lights_dir: Path,
        darks_dir: Optional[Path] = None,
        flats_dir: Optional[Path] = None,
        biases_dir: Optional[Path] = None,
        output_dir: Optional[Path] = None,
    ):
        self.lights_dir = lights_dir.resolve()
        self.darks_dir = darks_dir.resolve() if darks_dir else None
        self.flats_dir = flats_dir.resolve() if flats_dir else None
        self.biases_dir = biases_dir.resolve() if biases_dir else None
        self.working_dir = (output_dir or lights_dir.parent / "photometry").resolve()
        self.process_dir = self.working_dir / "process"
        self.masters_dir = self.working_dir / "masters"
        self.results_dir = self.working_dir / "results"

    # ── Initialisation ────────────────────────────────────────────────────────

    def _mkdirs(self) -> None:
        for d in [self.process_dir, self.masters_dir, self.results_dir]:
            d.mkdir(parents=True, exist_ok=True)

    def _count_fits(self, directory: Path) -> int:
        return len([f for f in directory.iterdir()
                    if f.suffix.lower() in (".fit", ".fits", ".fts")])

    # ── Étape 1 : Calibration et alignement ──────────────────────────────────

    def preprocess(self) -> str:
        """
        Exécute la pipeline de prétraitement Siril :
          link → (calibrate) → seqplatesolve → seqapplyreg

        Retourne le préfixe de la séquence enregistrée (ex: "r_pp_light_").
        """
        self._mkdirs()

        n = self._count_fits(self.lights_dir)
        if n == 0:
            raise FileNotFoundError(f"Aucun fichier FITS dans {self.lights_dir}")
        log.info("%d frames FITS dans %s", n, self.lights_dir)

        cmds: List[str] = []
        focal = SEESTAR_S30_PRO["focal_mm"]
        pixsz = SEESTAR_S30_PRO["pixel_um"]

        # 1a. Lien symbolique des frames dans process/
        cmds += [
            f'cd "{self.lights_dir}"',
            f'link light -out="{self.process_dir}"',
            f'cd "{self.process_dir}"',
        ]

        seq = "light_"

        # 1b. Calibration (si frames disponibles)
        has_darks = self.darks_dir and self._count_fits(self.darks_dir) > 0
        has_flats = self.flats_dir and self._count_fits(self.flats_dir) > 0
        has_biases = self.biases_dir and self._count_fits(self.biases_dir) > 0

        if has_darks:
            cmds += [
                f'cd "{self.darks_dir}"',
                f'link dark -out="{self.masters_dir}"',
                f'cd "{self.masters_dir}"',
                "stack dark rej 3 3 -nonorm -out=darks_master",
                f'cd "{self.process_dir}"',
            ]
        if has_flats:
            cmds += [
                f'cd "{self.flats_dir}"',
                f'link flat -out="{self.masters_dir}"',
                f'cd "{self.masters_dir}"',
                "stack flat rej 3 3 -norm=mul -out=flats_master",
                f'cd "{self.process_dir}"',
            ]
        if has_biases:
            cmds += [
                f'cd "{self.biases_dir}"',
                f'link bias -out="{self.masters_dir}"',
                f'cd "{self.masters_dir}"',
                "stack bias rej 3 3 -nonorm -out=biases_master",
                f'cd "{self.process_dir}"',
            ]

        if has_darks or has_flats or has_biases:
            cal = f"calibrate {seq} -cfa -equalize_cfa"
            if has_darks:
                cal += f' -dark="{self.masters_dir}/darks_master"'
            if has_flats:
                cal += f' -flat="{self.masters_dir}/flats_master"'
            if has_biases:
                cal += f' -bias="{self.masters_dir}/biases_master"'
            cmds.append(cal)
            seq = "pp_light_"
            log.info("Calibration activée → séquence : %s", seq)
        else:
            log.info("Pas de frames de calibration — traitement direct")

        # 1c. Résolution astrométrique de chaque frame
        cmds.append(
            f"seqplatesolve {seq} -nocache -force "
            f"-focal={focal} -pixelsize={pixsz} "
            f"-radius={SEESTAR_S30_PRO['fov_deg']:.1f}"
        )

        # 1d. Alignement par WCS
        cmds.append(
            f"seqapplyreg {seq} -framing=max -filter-round=2.5k"
        )

        registered = f"r_{seq}"

        if not run_siril_script(self.process_dir, cmds, "_preprocess.ssf"):
            raise RuntimeError(
                "Échec du prétraitement — consultez les logs Siril ci-dessus.\n"
                "Assurez-vous que siril-cli 1.4+ est installé et que les frames sont valides."
            )

        log.info("Prétraitement terminé. Séquence alignée : %s", registered)
        return registered

    # ── Étape 2 : Photométrie ─────────────────────────────────────────────────

    def run_photometry(
        self,
        seq_prefix: str,
        target_ra: float,
        target_dec: float,
        ref_stars: List[Tuple[float, float]],
        channel: int = PHOT_DEFAULTS["channel"],
    ) -> Path:
        """
        Configure la photométrie et génère la courbe de lumière.

        Paramètres :
            seq_prefix  Préfixe de la séquence alignée (ex: "r_pp_light_")
            target_ra   RA de l'étoile cible en degrés décimaux
            target_dec  Dec de l'étoile cible en degrés décimaux
            ref_stars   Liste de (RA, Dec) des étoiles de comparaison (min 3)
            channel     0=R  1=G (défaut)  2=B

        Retourne le chemin vers light_curve.dat
        """
        if len(ref_stars) < 3:
            log.warning(
                "%d étoiles de référence fournies — recommandé : 4+. "
                "La précision photométrique sera limitée.",
                len(ref_stars)
            )

        phot = PHOT_DEFAULTS
        cmds = [
            f'cd "{self.process_dir}"',
            # Paramètres de l'ouverture photométrique
            f"setphot "
            f"-aperture={phot['aperture_px']} "
            f"-inner={phot['inner_px']} "
            f"-outer={phot['outer_px']} "
            f"-gain={SEESTAR_S30_PRO['gain_e_adu']}",
        ]

        # Commande light_curve avec coordonnées WCS
        # -autoring : ajuste automatiquement les anneaux selon la FWHM mesurée
        lc_cmd = (
            f"light_curve {seq_prefix} {channel} -autoring "
            f"-wcs={wcs_arg(target_ra, target_dec)}"
        )
        for ref_ra, ref_dec in ref_stars:
            lc_cmd += f" -refwcs={wcs_arg(ref_ra, ref_dec)}"
        cmds.append(lc_cmd)

        if not run_siril_script(self.process_dir, cmds, "_photometry.ssf"):
            raise RuntimeError(
                "Échec de la photométrie.\n"
                "Vérifiez que :\n"
                "  - Les coordonnées de l'étoile cible sont dans le champ\n"
                "  - Les frames sont bien plate-solvées (WCS présent)\n"
                "  - Siril 1.4+ est installé"
            )

        # Siril écrit light_curve.dat dans le répertoire courant (process_dir)
        lc_src = self.process_dir / "light_curve.dat"
        if not lc_src.exists():
            raise FileNotFoundError(
                "light_curve.dat introuvable après l'exécution de Siril.\n"
                "L'étoile cible était peut-être hors champ ou non détectée."
            )

        lc_dest = self.results_dir / "light_curve.dat"
        shutil.copy2(lc_src, lc_dest)
        log.info("Courbe de lumière : %s", lc_dest)

        png_src = self.process_dir / "light_curve.png"
        if png_src.exists():
            png_dest = self.results_dir / "light_curve.png"
            shutil.copy2(png_src, png_dest)
            log.info("Graphe : %s", png_dest)

        return lc_dest

    # ── Étape 3 : Export AAVSO ────────────────────────────────────────────────

    def export_aavso(
        self,
        lc_dat: Path,
        target_name: str,
        observer_code: str = "XXXX",
        filter_name: str = "TG",  # TG = Transformed Green (standard pour OSC)
    ) -> Path:
        """
        Convertit light_curve.dat en CSV format AAVSO Extended.

        Note : les magnitudes sont DIFFÉRENTIELLES (relative aux étoiles de
        comparaison). Pour une soumission AAVSO complète, calibrez en magnitude
        absolue en utilisant les magnitudes V connues de vos étoiles de référence.

        Format light_curve.dat attendu (colonnes séparées par espaces/tabs) :
            #frame_index  JD  diff_mag  mag_err  [SNR]  [...]
        """
        rows: List[Tuple[float, float, float, float]] = []

        with open(lc_dat, encoding="utf-8") as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                try:
                    # Siril peut mettre l'index de frame ou la date JD en col 0
                    # On tente de déterminer la colonne JD (valeur > 2400000)
                    if float(parts[0]) > 2_400_000:
                        jd, mag, err = float(parts[0]), float(parts[1]), float(parts[2])
                        snr = float(parts[3]) if len(parts) > 3 else 0.0
                    else:
                        jd, mag, err = float(parts[1]), float(parts[2]), float(parts[3])
                        snr = float(parts[4]) if len(parts) > 4 else 0.0
                    rows.append((jd, mag, err, snr))
                except (ValueError, IndexError):
                    continue

        if not rows:
            log.warning("Aucune mesure valide dans %s", lc_dat)
            return lc_dat

        safe_name = target_name.replace(" ", "_").replace("/", "-")
        out_path = self.results_dir / f"{safe_name}_aavso.csv"

        with open(out_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            # En-têtes de directives AAVSO
            w.writerow(["#TYPE=EXTENDED"])
            w.writerow([f"#OBSCODE={observer_code}"])
            w.writerow(["#SOFTWARE=Siril + seestar_variable_photometry.py"])
            w.writerow(["#DELIM=,"])
            w.writerow(["#DATE=JD"])
            w.writerow(["#OBSTYPE=CCD"])
            # Colonnes
            w.writerow([
                "NAME", "DATE", "MAG", "MERR", "FILT", "TRANS",
                "MTYPE", "CNAME", "CMAG", "KNAME", "KMAG",
                "AMASS", "GROUP", "CHART", "NOTES",
            ])
            for jd, mag, err, snr in rows:
                w.writerow([
                    target_name, f"{jd:.6f}", f"{mag:.4f}", f"{err:.4f}",
                    filter_name, "NO", "STD",
                    "ENSEMBLE", "na",
                    "na", "na",
                    "na", "1", "na",
                    f"SNR={snr:.1f}|seestar_s30pro",
                ])

        log.info("Export AAVSO : %s  (%d observations)", out_path, len(rows))
        return out_path


# ─────────────────────────────────────────────────────────────────────────────
# Mode interactif
# ─────────────────────────────────────────────────────────────────────────────

def _ask_coords(label: str) -> Tuple[float, float]:
    """Demande interactivement les coordonnées d'une étoile."""
    print(f"\n  {label}")
    print("  Formats acceptés :")
    print("    Degrés décimaux : 83.8221,-5.3911")
    print("    Sexagésimal     : 05:35:17.3,-05:23:28")
    if HAS_SIMBAD:
        print("    Nom d'étoile    : RR Lyr, T Cep, eta Aql ...")
    while True:
        val = input("  > ").strip()
        if not val:
            continue
        if "," in val:
            try:
                return parse_radec(val)
            except Exception as exc:
                print(f"  [!] Erreur de format : {exc}")
        elif HAS_SIMBAD:
            result = resolve_star_name(val)
            if result:
                ra, dec = result
                print(f"  → RA={ra:.6f}°  Dec={dec:+.6f}°")
                return ra, dec
            print("  [!] Étoile introuvable dans Simbad.")
        else:
            print("  [!] Installez astroquery pour la résolution par nom. pip install astroquery")


def interactive_mode(lights_dir: Path) -> Dict:
    """Guide l'utilisateur pour configurer une session de photométrie."""
    print()
    print("=" * 65)
    print("   SEESTAR S30 PRO  —  Photométrie d'étoiles variables")
    print("=" * 65)
    print(f"\n   Répertoire lights : {lights_dir}")
    print(
        "\n   CONSEIL : Utilisez l'outil AAVSO Variable Star Plotter (VSP)\n"
        "   pour identifier votre étoile cible et ses étoiles de comparaison.\n"
        "   → https://www.aavso.org/apps/vsp/"
    )

    # Étoile cible
    print("\n── Étoile cible ──────────────────────────────────────────────")
    name = input("  Nom (ex: 'T Cep', 'RR Lyr', 'eta Aql') : ").strip()
    if not name:
        name = "target"
    ra, dec = _ask_coords(f"Coordonnées de {name!r}")

    # Étoiles de comparaison
    print("\n── Étoiles de comparaison ────────────────────────────────────")
    print("  Choisissez 3 à 6 étoiles stables dans le même champ.")
    print("  Critères AAVSO : non-variables, ±1 mag de la cible,")
    print("  couleur similaire, dans le même champ de vue (≈1.5°).\n")

    refs: List[Tuple[float, float]] = []
    i = 1
    while True:
        print(f"  Étoile de comparaison #{i} (laisser vide pour terminer) :")
        try:
            val = input("  > ").strip()
        except (KeyboardInterrupt, EOFError):
            break
        if not val:
            if len(refs) >= 3:
                break
            print("  [!] Minimum 3 étoiles de comparaison requises.")
            continue
        # Si l'utilisateur entre juste un nom
        if "," not in val and HAS_SIMBAD:
            result = resolve_star_name(val)
            if result:
                refs.append(result)
                print(f"  → RA={result[0]:.6f}°  Dec={result[1]:+.6f}°")
                i += 1
                continue
        try:
            ra_r, dec_r = parse_radec(val)
            refs.append((ra_r, dec_r))
            i += 1
        except Exception as exc:
            print(f"  [!] Format incorrect : {exc}")

    obs_code = input("\n  Code observateur AAVSO (laisser vide = XXXX) : ").strip() or "XXXX"

    return {
        "target_name": name,
        "target_ra": ra,
        "target_dec": dec,
        "ref_stars": refs,
        "observer_code": obs_code,
    }


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Courbes de lumière pour étoiles variables — Seestar S30 Pro + Siril",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemples :
  # Mode interactif (recommandé pour démarrer)
  python seestar_variable_photometry.py --lights /data/seestar/lights/

  # Avec étoile cible par nom + étoiles de comparaison en coordonnées
  python seestar_variable_photometry.py \\
    --lights /data/lights/ \\
    --target "RR Lyr" \\
    --refs "19:25:27.9,+42:47:04" "19:26:10.2,+42:50:45" "19:24:33.1,+43:01:23"

  # Reprendre la photométrie sans refaire le prétraitement
  python seestar_variable_photometry.py \\
    --lights /data/lights/ \\
    --config rr_lyr.json \\
    --skip-preprocess --seq-prefix r_pp_light_

  # Sauvegarder la configuration pour les prochaines sessions
  python seestar_variable_photometry.py \\
    --lights /data/lights/ \\
    --target "eta Aql" \\
    --save-config eta_aql.json
""",
    )
    p.add_argument("--lights", required=True,
                   help="Répertoire des frames FITS du Seestar (lights/)")
    p.add_argument("--darks", help="Répertoire des darks (optionnel)")
    p.add_argument("--flats", help="Répertoire des flats (optionnel)")
    p.add_argument("--biases", help="Répertoire des biases (optionnel)")
    p.add_argument("--output", help="Répertoire de sortie (défaut : lights/../photometry/)")
    p.add_argument("--target",
                   help="Étoile cible : nom Simbad (ex: 'RR Lyr') ou 'RA,Dec'")
    p.add_argument("--refs", nargs="+", metavar="RA,Dec",
                   help="Étoiles de comparaison 'RA,Dec' — répéter pour chaque étoile (min 3)")
    p.add_argument("--config", help="Charger la config étoiles depuis un fichier JSON")
    p.add_argument("--save-config", metavar="FILE",
                   help="Sauvegarder la config étoiles dans un fichier JSON")
    p.add_argument("--channel", type=int, default=PHOT_DEFAULTS["channel"],
                   choices=[0, 1, 2],
                   help="Canal couleur : 0=R  1=G (défaut)  2=B")
    p.add_argument("--observer", default="XXXX",
                   help="Code observateur AAVSO (ex: JXXX)")
    p.add_argument("--filter", default="TG",
                   help="Filtre AAVSO (défaut: TG = Transformed Green)")
    p.add_argument("--skip-preprocess", action="store_true",
                   help="Sauter le prétraitement (séquence déjà alignée)")
    p.add_argument("--seq-prefix", default=None,
                   help="Préfixe de séquence à utiliser avec --skip-preprocess")
    p.add_argument("--debug", action="store_true", help="Activer les logs de débogage")
    return p


def main() -> None:
    args = build_parser().parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    # ── Chemins ──────────────────────────────────────────────────────────────
    lights_dir = Path(args.lights)
    if not lights_dir.is_dir():
        log.error("Répertoire lights introuvable : %s", lights_dir)
        sys.exit(1)

    darks_dir = Path(args.darks) if args.darks else None
    flats_dir = Path(args.flats) if args.flats else None
    biases_dir = Path(args.biases) if args.biases else None
    output_dir = Path(args.output) if args.output else None

    # ── Configuration des étoiles ─────────────────────────────────────────────
    star_cfg: Optional[Dict] = None

    if args.config:
        cfg_path = Path(args.config)
        if cfg_path.exists():
            with open(cfg_path, encoding="utf-8") as fh:
                star_cfg = json.load(fh)
            # Convertir les refs de liste de listes en liste de tuples
            if "ref_stars" in star_cfg:
                star_cfg["ref_stars"] = [tuple(r) for r in star_cfg["ref_stars"]]
            log.info("Configuration chargée depuis %s", cfg_path)

    if star_cfg is None:
        if args.target:
            if "," in args.target:
                t_ra, t_dec = parse_radec(args.target)
                t_name = args.target
            else:
                t_name = args.target
                coords = resolve_star_name(t_name)
                if not coords:
                    log.error("Impossible de résoudre '%s' via Simbad.", t_name)
                    log.error("Entrez les coordonnées directement avec --target 'RA,Dec'")
                    sys.exit(1)
                t_ra, t_dec = coords

            if not args.refs or len(args.refs) < 3:
                log.error(
                    "Fournissez au moins 3 étoiles de comparaison avec --refs 'RA,Dec' ..."
                )
                sys.exit(1)

            refs = [parse_radec(r) for r in args.refs]
            star_cfg = {
                "target_name": t_name,
                "target_ra": t_ra,
                "target_dec": t_dec,
                "ref_stars": refs,
                "observer_code": args.observer,
            }
        else:
            # Mode interactif
            star_cfg = interactive_mode(lights_dir)

    # Sauvegarde optionnelle de la config
    if args.save_config:
        save_path = Path(args.save_config)
        with open(save_path, "w", encoding="utf-8") as fh:
            json.dump(star_cfg, fh, indent=2)
        log.info("Configuration sauvegardée dans %s", save_path)

    # ── Initialisation du pipeline ────────────────────────────────────────────
    pipeline = PhotometryPipeline(
        lights_dir=lights_dir,
        darks_dir=darks_dir,
        flats_dir=flats_dir,
        biases_dir=biases_dir,
        output_dir=output_dir,
    )

    # ── Prétraitement ─────────────────────────────────────────────────────────
    if args.skip_preprocess:
        seq_prefix = args.seq_prefix or "r_light_"
        log.info("Prétraitement ignoré — séquence utilisée : %s", seq_prefix)
    else:
        print(
            "\n[1/3] Prétraitement des frames"
            " (calibration + plate solve + alignement) ..."
        )
        seq_prefix = pipeline.preprocess()

    # ── Photométrie ───────────────────────────────────────────────────────────
    print("\n[2/3] Photométrie d'ouverture + courbe de lumière ...")
    lc_path = pipeline.run_photometry(
        seq_prefix=seq_prefix,
        target_ra=star_cfg["target_ra"],
        target_dec=star_cfg["target_dec"],
        ref_stars=star_cfg["ref_stars"],
        channel=args.channel,
    )

    # ── Export AAVSO ──────────────────────────────────────────────────────────
    print("\n[3/3] Export AAVSO ...")
    aavso_path = pipeline.export_aavso(
        lc_dat=lc_path,
        target_name=star_cfg["target_name"],
        observer_code=star_cfg.get("observer_code", args.observer),
        filter_name=args.filter,
    )

    # ── Résumé final ──────────────────────────────────────────────────────────
    png_path = pipeline.results_dir / "light_curve.png"
    print()
    print("=" * 65)
    print("   TERMINÉ")
    print("=" * 65)
    print(f"   Étoile cible       : {star_cfg['target_name']}")
    print(f"   RA / Dec           : {star_cfg['target_ra']:.4f}° / {star_cfg['target_dec']:+.4f}°")
    print(f"   Étoiles de réf.    : {len(star_cfg['ref_stars'])}")
    print(f"   Canal photométrie  : {['R','G','B'][args.channel]} (ch {args.channel})")
    print()
    print(f"   Données brutes     : {lc_path}")
    if png_path.exists():
        print(f"   Graphe             : {png_path}")
    print(f"   Export AAVSO       : {aavso_path}")
    print()
    print("   Prochaines étapes :")
    print("   1. Vérifiez la courbe dans le graphe PNG")
    print("   2. Calibrez en magnitude absolue si nécessaire")
    print("      (ajoutez les magnitudes V connues de vos étoiles de comparaison)")
    print("   3. Soumettez sur https://www.aavso.org/webobs")
    print("=" * 65)


if __name__ == "__main__":
    main()
