# SageMath algebraic-geometry audit for the one-loop cusp constraints.
# Homogeneous variables for line Plücker coordinates can be substituted as needed.
R.<p12,p23,p34,p14,p13,p24> = PolynomialRing(QQ)
I = ideal([p12, p23])
print("Cusp ideal:", I)
print("Dimension:", I.dimension())
print("Radical:", I.radical())
print("Projective closure / primary decomposition should be run in the chosen Plücker chart.")
