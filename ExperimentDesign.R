my_data=read.csv("https://raw.githubusercontent.com/Harrindy/StatEngine/master/Data/TensileStrength.csv")
head(my_data,2)
my_data$Concentration=as.factor(my_data$Concentration)
boxplot(Strength~Concentration,data=my_data)
res.aov=aov(Strength~Concentration,data=my_data)
summary(res.aov)
TukeyHSD(res.aov)
TukeyHSD(res.aov,conf.level=0.95)

my_data=read.csv("https://raw.githubusercontent.com/Harrindy/StatEngine/master/Data/AdhesionForce.csv")
head(my_data,2)
my_data$Primer=as.factor(my_data$Primer)
my_data$Method=as.factor(my_data$Method)
head(my_data,2)
boxplot(Adhesion~Primer,data=my_data)
boxplot(Adhesion~Method,data=my_data)
boxplot(Adhesion~Primer:Method,data=my_data)
res.aov=aov(Adhesion~Primer+Method+Primer:Method,data=my_data)
summary(res.aov)
TukeyHSD(res.aov,"Primer")
TukeyHSD(res.aov,"Method")
TukeyHSD(res.aov,"Primer:Method")



my_data=read.csv("https://raw.githubusercontent.com/Harrindy/StatEngine/master/Data/CodedSurfaceRoughness.csv")

# We need to tell R all the factors
head(my_data,2)
my_data$Feed=as.factor(my_data$Feed)
my_data$Depth=as.factor(my_data$Depth)
my_data$Angle=as.factor(my_data$Angle)
boxplot(Roughness~Feed,data=my_data)
boxplot(Roughness~Depth,data=my_data)
boxplot(Roughness~Angle,data=my_data)

# Run ANOVA to find P-value fot the tests
res.aov=aov(Roughness~Feed+Depth+Angle+
                Feed:Depth+Feed:Angle+Depth:Angle+
                Feed:Depth:Angle,data=my_data)
summary(res.aov)

# If rejecting H0, TukeyHSD finds which pairs caused the difference.
TukeyHSD(res.aov,'Feed')
