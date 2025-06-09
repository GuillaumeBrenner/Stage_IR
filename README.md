# zibell

**zibell** est un package Python qui implémente la simulation et l’estimation robuste d’un modèle **Zero-Inflated Bell (ZI-Bell)** en utilisant la méthode d’estimation **MDPDE (Minimum Density Power Divergence Estimator)**.

---

## 📘 Description

Ce package permet de :

- Générer des données simulées selon un modèle ZI-Bell ;
- Estimer les paramètres du modèle par la méthode MDPDE ;
- Réaliser des expériences de Monte Carlo pour évaluer les performances de l’estimation.

Le modèle ZI-Bell est un modèle à inflation de zéros où les données sont composées d’une proportion non nulle de zéros, et les valeurs non nulles suivent une loi de Bell. Le paramètre d'inflation est modélisé par une fonction logistique.

---

## 📦 Installation

Assurez-vous d’avoir Python ≥ 3.8.

Clonez le dépôt, puis installez en mode développement :

```bash
cd zibellpackage
pip install -e .
```

---

## 🚀 Utilisation

Voici un exemple simple d'utilisation :

```python
from zibell.simulate import data
from zibell.estimation import estime_rep

# Génération de données
n = 100
b = [-0.5, 1.2]  # Coefficients de la moyenne (mu)
g = -1.1         # Coefficient d'inflation des zéros
alpha = 0.1      # Paramètre de robustesse

donnees = data(n, b, g)
Y, X, Z = donnees['Y'], donnees['X'], donnees['Z']

# Estimation des paramètres
theta = estime_rep(Y, X, Z, alpha)
print("Paramètres estimés :", theta)
```

---

## 🧠 Contenu du package

- `simulate.py` : génération des données ZI-Bell ;
- `estimation.py` : estimation des paramètres via MDPDE ;
- `metrics.py` *(optionnel)* : outils pour calculer biais, RMSE, etc.

---

## 📊 Méthodologie

L’estimation repose sur la fonction de divergence de puissance (Power Divergence) entre les densités théoriques et empiriques. Cela permet une robustesse accrue aux données aberrantes (outliers).

---

## ✅ Dépendances

- `numpy`
- `scipy`
- `lambertw` (si utilisé pour l’inversion de `W`)
- `tqdm` (optionnel pour les boucles Monte Carlo)

---

## 📁 Structure recommandée

```
zibell/
│
├── zibell/
│   ├── __init__.py
│   ├── simulate.py
│   ├── estimation.py
│   └── metrics.py (optionnel)
│
├── test_usage.py
├── setup.py
└── README.md
```

---

## ✍️ Auteur

Projet initialement développé en R puis porté en Python dans le cadre de l'Iniatiation à la recherche.
HOUNKPATI Codjovi Guillaume, Développeur en alternance
Dr. Solym Manou-Abi, Maitre de conférences, Laboratoire de Mathématiques et Applications
CNRS UMR 7348 – Université de Poitiers

---
