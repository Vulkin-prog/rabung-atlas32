# Publication successive des quatre articles

L'ordre retenu est : **nouvelle édition des certificats, méthode GPU, atlas,
quadruplets cubiques**. Le premier est une révision d'un article déjà publié ;
les trois suivants sont de nouveaux articles. La série conserve donc quatre
articles, avec plusieurs éditions possibles du premier.

Une publication désigne ici une édition scientifique rendue publiquement
accessible et citable, par dépôt ou par revue. Un dépôt public et une évaluation
par les pairs ont des statuts distincts. Les certificats disposent de deux
éditions publiées : v1.0.0, DOI **10.5281/zenodo.22980965**, et v2.0,
DOI **10.5281/zenodo.23337824**, publiée le 13 octobre 2026. La méthode GPU
v1.0 est publiée le même jour, DOI **10.5281/zenodo.23338707**, ainsi que
l'atlas v1.0, DOI **10.5281/zenodo.23339668**. Seul le papier des quadruplets
cubiques attend son dépôt ; aucun identifiant futur n'est supposé.

| Étape | Article | Justification et contenu |
| --- | --- | --- |
| 1 | Certificats et bornes, nouvelle édition | Les six nouveaux certificats et leurs quatorze conséquences directes utilisent des scans CPU complets. Leur validité ne dépend pas de la publication du scanner ou de l'atlas. Cette édition porte la preuve du critère, les vérifications ponctuelles, les produits à module exact et les comparaisons par récurrences. |
| 2 | Méthode GPU | Le théorème de couverture et de complétion, les garanties arithmétiques et les deux comparaisons contrôlées forment un sujet autonome. L'application aux constructions cite les éditions publiées des certificats, en distinguant les pièces originales de la preuve et des résultats actuels. |
| 3 | Atlas arithmétique | Il peut alors citer la méthode du scanner et les certificats. Il conserve le recensement, les densités de Kummer, les normes, les codes locaux et les diagnostics ternaires. Les résultats elliptiques cubiques détaillés sont réservés à l'article suivant. |
| 4 | Quadruplets cubiques | L'atlas fournit un contexte publié pour les questions sur les structures et les fluctuations. La formule exacte, les six courbes, les torsions CM et les quatre moments restent démontrés dans cet article. Le raccord analytique doit être contre-revu avant sa publication. |

## Références disponibles et ajouts par étape

Le papier des certificats v2.0 conserve sa référence à la v1.0.0. La méthode
GPU et l'atlas citent la v2.0 pour le critère et les résultats actuels, et
conservent la v1.0.0 pour ses pièces et comparaisons originales. Le papier
cubique cite la v2.0 pour le contexte des constructions finies. L'atlas cite
aussi la méthode GPU v1.0 pour la couverture, les hypothèses d'exactitude
et la distinction entre comparaisons contrôlées et temps du recensement.
Le papier cubique cite l'atlas v1.0 pour le contexte du recensement fini
et distingue ses diagnostics des moments obtenus en moyenne sur les premiers.
Les bibliographies ne citent aucun des articles encore non publiés. Les
références scientifiques extérieures sont conservées. Cette règle s'applique
aux références entre ces quatre articles ;
elle n'exclut pas les prépublications publiques extérieures déjà citées.

Les insertions prévues pour la première diffusion des quatre articles sont
effectuées : certificats v2.0 vers les trois articles suivants, méthode GPU
v1.0 vers l'atlas, puis atlas v1.0 vers les quadruplets cubiques.
Le dernier DOI reçu est celui de l'atlas, **10.5281/zenodo.23339668**.

Le papier cubique ne doit citer la méthode GPU que si une contribution précise
de celle-ci est effectivement utilisée. Ses preuves actuelles et ses certificats
elliptiques CPU ne l'exigent pas. La cohérence de la série repose sur des rôles
explicites, avec les mêmes conventions, plutôt que sur une bibliographie identique.

Les valeurs historiques restent attribuées à leur édition exacte. Par exemple,
la comparaison avec la borne publiée de 31 385 622 833 conserve sa référence
v1.0.0 ; la preuve d'une nouvelle borne peut recevoir une référence distincte
à la nouvelle édition. Une nouvelle édition ne justifie pas le remplacement
global de tous les anciens DOI.

## Contenu des textes et autonomie

Les titres scientifiques restent stables. Chaque introduction distingue son
objet : validité d'une construction, exactitude d'un calcul, structure d'une
population finie ou fluctuations d'un compte de longueur fixe. Tous les textes
emploient `W(r,k)`, couleurs d'abord, et distinguent `M`, `L` et `K`.
Les populations probabilistes et les domaines des calculs restent explicites.

L'atlas ne présente plus les constantes des quatre moments ni la formule
elliptique cubique d'un article futur. Son propre compte quadratique de
quadruplets et ses diagnostics ternaires restent présents. Le papier cubique
conserve l'intégralité de ses résultats et de leurs preuves. Une référence à un
article déjà publié peut situer une contribution sans remplacer les définitions,
corrections de bord et arguments nécessaires à la compréhension du nouveau texte.

## Procédure à chaque publication

1. Préparer un paquet public propre à l'article : son PDF, sa source, son code,
   ses données et certificats nécessaires, ses droits et sa notice de reproduction.
   Les autres brouillons, revues reçues et historiques éditoriaux restent dans
   le dossier de travail. Les composants communs publiés sont identifiés par
   leurs versions et empreintes.
2. Vérifier ce paquet dans un dossier extrait neuf. La révision des certificats
   a notamment besoin de son registre, de ses sources CPU, de ses neuf résultats,
   du programme de fermeture et des entrées publiées importées. L'atlas ajoute
   les données et contrôles propres au recensement. L'archive intégrale du
   dépôt de travail n'est pas le paquet public d'une seule étape.
3. Publier l'édition, puis enregistrer son titre exact, sa version, sa date et son
   DOI ou autre adresse pérenne effectivement délivrée.
4. Ajouter les citations pertinentes aux textes suivants, avec un rôle précis
   et des locateurs de théorème ou de section lorsque nécessaires. Compiler et
   vérifier leurs bibliographies après ces ajouts.
5. Conserver l'édition publiée telle quelle. Une évolution scientifique de ce
   texte donne lieu à une nouvelle édition identifiable.

`docs/PUBLICATION_SEQUENCE.json` enregistre l'ordre, les éditions
publiées et les références présentes. `scripts/verify_publication_references.py`
contrôle les titres, versions, locateurs et l'absence de références anticipées.
Ce contrôle ne publie aucun fichier et ne déduit pas la publication d'un commit
privé. Le registre et les citations seront actualisés à réception des véritables
identifiants des prochaines éditions.

Pour préparer les paquets séparés, le vérificateur intégré propose `--scope atlas`
et `--scope cubic` ; son comportement par défaut conserve le contrôle complet.
L'option `--catalogue` appartient au périmètre de l'atlas. Les deux commandes
de certificats indiquées dans le papier cubique fonctionnent directement sur
son propre paquet. Le manifeste de fichiers livré avec chaque publication
devra lui aussi être limité au paquet effectivement livré.

## État des éditions au 13 octobre 2026

Les certificats v2.0, la méthode GPU v1.0 et l'atlas v1.0 sont publiés.
Leurs PDF et archives déposés sont conservés tels quels ; leurs empreintes
correspondent aux fichiers précédemment livrés. Le papier cubique cite
désormais l'atlas publié pour le contexte des questions arithmétiques,
en distinguant les diagnostics finis des moments en moyenne sur les premiers.
Son PDF et son archive sont régénérés après cet ajout. La dernière publication
de la séquence est celle des quadruplets cubiques. Aucune insertion de
référence à un article futur ne reste prévue. La v1.0.0 des certificats reste
citée lorsqu'une comparaison porte sur ses résultats ou ses pièces originales.

Les quatre archives sont produites par `scripts/build_publication_archives.py`.
Chaque archive possède son README, les pièces de son article, un manifeste et
ses commandes de contrôle. Le jeu de données complet de l'atlas demeure un
composant séparé ; son accès doit accompagner une diffusion publique de l'atlas.
