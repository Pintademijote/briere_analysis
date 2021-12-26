m1=glm(S~extinction_rate, data=TIB, family = "quasipoisson")
summary(m1)

m2=lm(C~extinction_rate, data=TIB)
summary(m2)

m2=lm(indegree~extinction_rate, data=TIB)
summary(m2)

m2=lm(outdegree~extinction_rate, data=TIB)
summary(m2)

m2=lm(basal~extinction_rate, data=TIB)
summary(m2)

m2=lm(top~extinction_rate, data=TIB)
summary(m2)

m2=lm(intermediate~extinction_rate, data=TIB)
summary(m2)

m2=lm(omnivory~extinction_rate, data=TIB)
summary(m2)
