# Notes techniques — pipeline photométrie

---

## Versions

| Version | Changements |
|---------|-------------|
| v0.1.5 | chemin PRIMARY via `findcompstars` + `-ninastars` ; V_app ensembliste ; FWHM depuis `.seq` |
| v0.1.4 | marge de sécurité bord de champ portée à 200 px |
| v0.1.3 | filtrage in-frame des étoiles de comparaison FALLBACK A/B |
| v0.1.2 | vérification précoce des frames à l'étape 5 ; suppression de `annotate` Siril |
| v0.1.1 | filtrage des frames par WCS |
| v0.1.0 | première version |

---

## v0.1.5 — Comment ça marche

### 1. Chemin de session

L'utilisateur donne le chemin du dossier de session, qui doit contenir un sous-dossier `lights/` avec les frames FITS du Seestar. Au chargement, le script lit les headers FITS du premier fichier pour récupérer :

- `OBJCTRA` / `OBJCTDEC` ou `RA` / `DEC` → centre du champ (coordonnées équatoriales)
- `NAXIS1` / `NAXIS2` + échelle de plaque → taille du champ en degrés
- Nombre de fichiers `.fit` dans `lights/`

### 2. Requête du catalogue VSX

Une recherche en cône est envoyée à **VizieR** (`B/vsx/vsx`) autour du centre du champ. Le rayon de recherche est la moitié de la diagonale du champ de vue, arrondie à 0.1° supérieur — assez large pour inclure les étoiles proches du bord.

Filtre par défaut : magnitude maximale **14**. Toutes les étoiles variables VSX dans ce rayon et sous ce seuil sont retournées et affichées dans le tableau. L'utilisateur peut filtrer par nom ou type dans l'interface, et ajuster la magnitude max.

### 3. Sélection de la cible et des étoiles de comparaison

Après sélection d'une cible dans le tableau, le script appelle `findcompstars` (Siril) ou VizieR APASS DR9 selon la disponibilité. Voir la section [sélection du chemin photométrique](#sélection-du-chemin-photométrique) pour le détail.

### 4. Filtre et frames de calibration

**Filtre :** choix entre LP (défaut Seestar) ou sans filtre. Affecte uniquement le champ `FILT` dans l'export AAVSO (`CV` dans les deux cas — pas de bande photométrique standard).

**Calibration :** l'utilisateur peut fournir des dossiers de darks, flats et/ou bias. Si au moins un type est fourni, le pipeline crée les masters automatiquement à l'étape 0 :

```
link dark_ → stack dark_ rej 3 3 -nonorm → master_dark.fit
link flat_  → stack flat_  rej 3 3 -norm=mul → master_flat.fit
link bias_  → stack bias_  rej 3 3 -nonorm → master_bias.fit
calibrate light_ -dark=... -flat=... -bias=...  → pp_light_*.fit
```

Le stacking utilise le rejet sigma Winsorized 3σ, robuste aux cosmiques et trainées satellite. Sans calibration, la séquence active reste `light_`.

### 5. Exécution — étapes et algorithmes

**Étape 1+2 — Registration**

```
link light -out=process/
register light_ -2pass
```

`register -2pass` fait du pattern matching d'étoiles, entièrement offline. La transformation calculée par frame est une similarité : translation + rotation + scale uniforme, pas de correction de distorsion. Les résultats (matrices de transformation) sont écrits dans `light_.seq`.

**Étape 3 — Application des transformations**

```
seqapplyreg light_ -framing=max -filter-round=2.5k  →  r_light_*.fit
```

`-framing=max` : les frames de sortie utilisent l'union de tous les champs, donc sont plus grandes que l'entrée. Les étoiles proches du bord apparaissent dans moins de frames — d'où la marge d'exclusion de 200 px pour les étoiles de comp (v0.1.4).

`-filter-round=2.5k` : garde les 2500 meilleures frames par élongation stellaire. Sur une session Seestar de 200–500 frames ça ne rejette généralement rien, mais élimine les frames avec vibration ou trainée.

**Étape 4 — Plate solve par frame**

```
seqplatesolve r_light_ -nocache -force -focal=160 -pixelsize=2.9 -radius=2.5
```

Chaque frame reçoit sa propre solution astrométrique contre le catalogue Gaia DR3 (en ligne). Les headers WCS (`CRVAL`, `CRPIX`, matrice `CD`) sont écrits dans chaque `r_light_*.fit`. Nécessaire parce qu'après `-framing=max`, le WCS de la frame de référence n'est pas valide pour les autres.

**Étape 5 — Photométrie**

```
findcompstars <target> -catalog=APASS -narrowband=0 -max_stars=N  →  comp_stars.csv
setphot -aperture=10 -inner=20 -outer=30 -gain=1.0
light_curve r_light_ 0 -ninastars=comp_stars.csv
```

`light_curve` fait de la photométrie d'ouverture sur chaque frame. Le canal `0` correspond à la seule couche des frames mono issues de la registration CFA. La sortie brute est `light_curve.dat` (JD, V-C, erreur).

**Post-traitement Python**

Après la sortie de Siril, le pipeline Python :
1. Lit `light_curve.dat`
2. Calcule `V_app = V_C + median(V_catalog_comp)` depuis `comp_stars.csv`
3. Extrait le FWHM par frame depuis `r_light_.seq` (lignes `R0`), aligne avec `DATE-OBS` des FITS, convertit en arcsec
4. Écrit les fichiers de résultats

### 6. Sorties

```
results/StarName/
├── photometry.csv       JD, V_C, V_app, err
├── fwhm.csv             JD, FWHM_x_arcsec, FWHM_y_arcsec
└── StarName_aavso.csv   AAVSO Extended Format
```

`photometry.csv` :
```
# Ensemble V (APASS comp stars median): 12.284
# V_app = V_C + ensemble_V
JD,V_C,V_app,err
2461161.334643,1.2421,13.5261,0.0868
```

`fwhm.csv` :
```
# Plate scale: 1.0350 arcsec/px
JD,FWHM_x_arcsec,FWHM_y_arcsec
2461161.334643,4.520,4.310
```

---

## Sélection du chemin photométrique

**PRIMARY** — `findcompstars` + `-ninastars`

Siril gère lui-même la conversion coordonnées célestes → pixels via les WCS. Nécessite WCS valide sur toutes les frames et au moins 3 comp stars dans le champ.

**FALLBACK A** — coordonnées pixel depuis VizieR

Comp stars récupérées sur VizieR APASS DR9 (`e_Vmag < 0.05`, `|delta_V| < 2.0`), converties en pixel par projection TAN depuis le WCS de la frame de référence. Coordonnées arrondies à l'entier (Siril 1.4.x refuse les décimales).

**FALLBACK B** — une seule étoile proche

Dernier recours. Magnitude différentielle plus bruitée, variabilité intrinsèque de la comp star non détectable.

---

## V_app : biais filtre LP

`V_app = V_C + median(V_catalog)`

V_C n'est pas affecté par le filtre LP (cible et comp passent par le même filtre). Le biais entre en jeu via `median(V_catalog)` qui est en V standard APASS. Pour une nova naine (continuum bleu + compagnon rouge) vs des comp G/K, attendre ±0.2–0.5 mag en absolu. Les variations relatives dans une session sont fiables.

`seqsetmag` (calibration magnitude dans Siril) n'est pas scriptable en headless — V_app est donc calculé externement.

---

## Bugs Siril 1.4.x contournés

| Bug / limite | Contournement |
|-------------|---------------|
| Coords pixel en `-at` doivent être des entiers | `round()` avant de construire la commande |
| `-autoring` incompatible avec `-at` | `setphot` avec valeurs fixes |
| `seqsetmag` non scriptable en headless | V_app calculé depuis `comp_stars.csv` |
| Parser `.ssf` inclut les guillemets dans le path | Toujours `-out=path` sans guillemets |
| `seqpsf` headless : sortie console uniquement | Non utilisé ; FWHM lu depuis `.seq` |
| `light_curve -wcs` peut échouer sur frames sans WCS | PRIMARY utilise `-ninastars` ; FALLBACK A utilise `-at` |

---

## Échelle de plaque

```
206.265 * 2.9 µm / 160 mm = 3.74 arcsec/px   (pixel Bayer brut)
```

La valeur utilisée pour FWHM `.seq` → arcsec est **1.035 arcsec/px** — le Seestar empile en interne des sub-expositions de 10s avec un léger drizzle, ce qui modifie légèrement l'échelle effective.

> TODO : vérifier 1.035 empiriquement depuis la matrice CD d'une frame plate-solvée.

---

## Gain e⁻/ADU

`get_gain_eadu()` cherche dans l'ordre : `EGAIN`, `EPERDN`, `GAIN_E`, `CCDGAIN`, puis `GAIN`. Valeur acceptée seulement si `0.05 < g < 30`. Le Seestar écrit `GAIN=200` (réglage caméra, pas e⁻/ADU) — rejeté par cette plage, retombe sur **1.0 e⁻/ADU**.

Valeur physique IMX585 au gain 200 : ~0.5–0.7 e⁻/ADU. Les barres d'erreur sont donc ~20–40% trop larges — sans impact sur la photométrie différentielle. L'app de pilotage Seestar devra écrire `EGAIN` dans les headers pour corriger ça.
