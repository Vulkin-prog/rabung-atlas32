# Réponse aux deux rapports reçus le 8 octobre 2026

La révision initiale après les rapports a produit **l'atlas v1.3 (29 pages)**
et le **compagnon cubique v1.1 (16 pages)**. Elle part de la branche de la PR #1 au commit
`59a9f9c501ab6781996222470049940ddef2e691` et conserve les résultats scientifiques
scellés. Les deux rapports ont été lus intégralement; les antécédents décisifs
ont été contrôlés dans les sources primaires. Le compte assemblé et ses
moments restent le sujet du compagnon. L'identité générale de corrélation
et la famille géométrique sont maintenant attribuées explicitement.

La correction éditoriale suivante produit **l'atlas v1.4 (29 pages)** et le
**compagnon v1.2 (16 pages)**. Conformément à la demande de l'auteur, les
reviews et documents internes sans destination de publication sont retirés
des sources scientifiques des articles. Les trois citations de la review
privée du 4 octobre, sa notice bibliographique, la notice des notes de projet
et celle du manuscrit GPU non fourni sont supprimées. Les renvois aux dossiers
de revue et les récits de révision quittent également le texte des articles.
Les preuves, identités et certificats sont présentés directement. Les deux
articles destinés à la publication conservent leurs références mutuelles;
les codes et données sont identifiés dans leurs sections de reproductibilité.
Ce document reste un registre de travail et n'est cité dans aucun manuscrit.

## Rapports d'entrée

| Document fourni | SHA256 |
| --- | --- |
| `AUDIT_PRIORITE_ATLAS_CUBIQUE_2026-10-08_FR.pdf` | `1990b4dceb74ae31e3aeca6e38231c9f2459b3ffc6de16337f574596202c08f1` |
| `CONTRE_REVUE_ATLAS_ET_QUADRUPLETS_2026-10-08_FR.pdf` | `bf4a63d595b8f4e24858eb69ad36faed31ab03d1b363a2b4c3c8f1b079e85fbd` |

Ces empreintes identifient les retours utilisés. Les PDF externes et les
archives privées de travail ne sont pas redistribués dans le dépôt.

## Corrections et décisions

| Retour | Modification apportée | Contrôle ou limite |
| --- | --- | --- |
| L'identité elliptique a un antécédent exact Greene–Lennon. | Résumé et introduction recentrés sur le compte et les moments. L'ancien théorème de corrélation devient un corollaire attribué. La section 2.2 donne le changement fractionnel, la définition normalisée, la phase, le point à l'infini et l'application du théorème 1.1 de Lennon. | Les sources Greene (définition 3.5 et théorème 4.4(i)) et Lennon ont été relues. Un contrôle indépendant de l'exposition couvre 2 360 quadruplets, dont A=0 et des phases non triviales. |
| La famille de genre deux, ses facteurs et l'isogénie sont déjà connus. | Freeman–Satoh, proposition 4.2, et Guillevic–Vergnaud, §3.1, p. 241, sont cités au point d'usage. La preuve géométrique directe demeure pour rendre explicites les twists et les contributions de bord. | Les modèles courts et la portée de la proposition ont été contrôlés dans les textes primaires. Aucune nouveauté générale d'isogénie n'est revendiquée. |
| Le passage du modèle fractionnaire au modèle entier est écrit dans le mauvais sens. | La section 3.2 du compagnon écrit la flèche `(X,Y) -> (9X,27Y)` de E(10/3) vers E3 et indique le coefficient constant 718/9 du modèle source. | Le transport du certificat dans l'oracle était déjà correct; les comptes et traces conservés ne changent pas. |
| L'argument sur les deux bords ordinaires est faux. | La preuve de propagation de l'atlas utilise seulement l'impossibilité de `w+j=0` pour j=0..31, puis en déduit que l'intervalle à droite du représentant w ne traverse pas p. Elle n'affirme plus w>=32. | Les identités et conclusions de propagation sont conservées; la majoration non nulle de Q(-j) suffit. |
| Les 100 modèles résiduels n'ont pas de qualification de corps rejouée. | Ajout de `RESIDUAL_IRREDUCIBILITY.json` et du vérificateur `residual_irreducibility.py`: 100 liaisons polynomiales, 100 témoins sur Q et 100 sur Q(i). `verify_synthesis.py` les exige avant le calcul des images. | Réductions de degré conservé, critère exact de Frobenius et pgcd. Les témoins sur Q(i) utilisent q=1 mod 4. Les plus grands q sont 59 et 181. Tests négatifs: réduction réductible, perte de degré, premier non scindé, polynôme différent et témoin altéré. |
| Les contrastes ternaires et leur centrage sont ambigus. | L'atlas définit les deux demi-millions, delta N, Delta, les innovations D(N) et D(T), puis les deux corrections iid epsilon(N) et epsilon(T). Le tableau donne les comptes bruts aux seuils 21, 24 et 27. | Les valeurs 217, 5296/27 et 563/27 sont distinguées des contrastes entièrement centrés. Les pourcentages de retrait précédant le centrage du filtre sont nommés ainsi. Le terme terminal des mots est borné, pas annulé par les seules tailles égales. |
| « Vertical » prête à confusion pour des courbes fixes et des premiers variables. | Les manuscrits et la documentation courante utilisent « moments over primes » ou « prime-averaged moments ». | Les noms et textes des archives scientifiques d'origine gardent leur identité documentaire. |
| Le raccord analytique doit être plus explicite. | Numéros exacts de Ramakrishnan et versions arXiv de Gillespie–Ji et Chen–Shen précisés. Regroupement (2,1,1) attribué à Ramakrishnan. Cartes des trois twists CM sur K(cuberoot(3)) et blocs de Weil–Deligne aux places ramifiées développés. Fité–Sutherland cité pour le conditionnement CM par classe finie. | Les conclusions restent les quatre moments du compte sur les véritables premiers. L'énumération des routes n'est pas une preuve automatique des théorèmes analytiques. |
| L'atlas doit distinguer théorie classique et réalisations propres. | Introduction centrée sur le recensement, les certificats Kummer spécifiques, la génération complète des relations par normes dans deux exemples, et le code arithmétique avec sa loi locale de suites. Wright, Karthick Babu et coauteurs, Delsarte et Cameron–Kusuma–Solé ajoutés. Kiritchenko et coauteurs cités à la quartique de longueur quatre. | Le facteur 4096 et le lien distance duale/indépendance limitée sont présentés comme des conséquences classiques du code certifié. Aucun « premier résultat mondial » n'est annoncé. |
| Références et détails éditoriaux à corriger. | Publications de Monroe, Gillespie–Ji et Chen–Shen ajoutées; titres publiés et titres des versions lues distingués. Monroe v7 reste la source historique explicitement comparée. Milne et PARI sont figés par version. « radicands », la ponctuation de K, le domaine Q[X] et « power-class group » sont corrigés. Les paramètres d'exécution passent en annexe. | La référence PARI renvoie à l'archive officielle 2.19.0 et à son empreinte, pas au manuel mobile html-stable. Les antécédents historiques des triples et quadruplets cubiques sont ajoutés au compagnon. |
| La revendication de calcul doit avoir une provenance publique immuable. | Ajout de `PUBLICATION_REPRODUCIBILITY.md`: carte des manuscrits, scanner gelé, données, certificats et commandes; distinction entre identifiants de versions et divulgation publique. | **Préparation réalisée, publication non réalisée.** Le dépôt reste privé et la PR en brouillon. Une archive publique persistante et les décisions de droits restent à organiser. Un hash privé ne prouve pas une date de divulgation publique. |

## Vérifications de la révision

- Les 15 tests unitaires passent, dont les quatre nouveaux tests d'irréductibilité
  et de liaison des polynômes.
- Les 200 témoins d'irréductibilité ont été générés puis rejoués. Par degré 2..6,
  les premiers maximaux sont `(41,29,31,23,59)` sur Q et
  `(61,89,61,61,181)` sur Q(i).
- La déduction hypergéométrique passe sur 2 360 cas jusqu'au premier 97:
  49 paramètres A=0 et 1 589 phases cubiques non triviales. Les sommes
  littérales, le transport fractionnel, les traces directes, les paramètres
  opposés et le modèle à 3-torsion sont comparés en arithmétique exacte.
- Le vérificateur algébrique complet passe avec les nouvelles qualifications:
  10 certificats multiracines, 445 identités de puissances quatrièmes,
  808 064 valeurs de caractères, 35 identités de normes, les 100 images
  résiduelles maximales et les lois locales existantes.
- Les deux sources compilent deux fois sans boîte débordante ni référence
  indéfinie. Les 45 pages ont été rendues et examinées; les pages contenant
  les nouvelles formules et tables ont aussi été contrôlées en pleine taille.
- Les reproductions éditoriales passent: 2 625 certificats elliptiques,
  207 petits premiers pour les moments, 315 routes analytiques finies et les
  contrôles exacts du même filtre. La reclassification CPU des 1 234 808
  positions ordinaires reproduit les deux listes archivées à l'octet près.
- Le contrôle complémentaire `verify_report4.py` passe. Les 48 identités de
  sources copiées sont conservées et les 70 empreintes éditoriales concordent.

La CI exige ces contrôles sur la révision commise, ainsi que les contrôles
historiques du dépôt et la compilation des manuscrits. Le manifeste éditorial
est renouvelé séparément des sceaux scientifiques originaux. La vérification
locale ne relance pas le recensement GPU complet.

## Portée restante

Les constantes des moments et les comptes exacts ne sont pas changés par
cette révision. Les deux retours aident à contrôler les preuves et la priorité;
ils ne constituent ni une revue humaine spécialisée de journal ni une
formalisation dans un assistant de preuve. Une contre-relecture humaine du
raccord analytique reste recommandée avant soumission. Les quatre moments
ne prouvent ni la loi limite complète, ni un centrage brut à erreur o(1),
ni l'intensité critique ou la probabilité de vide à longueur logarithmique.
