import math

# get info from problem
I = int(input("Number of levels: "))
print(I)

J = int(input("Sample size of the "+str(I)+" levels: "))
print(J)

print("List the", I, "sample means: ")
ybari = []
while(len(ybari) < I):
    ybari.append(float(input("")))
print(ybari)

print("List the", I, "sample variances: ")
Syi = []
while(len(Syi) < I):
    Syi.append(float(input("")))
print(Syi)

# compute ANOVA values
anova = [["Src", "dof", "SS", "MS", "F"],
         ["Tmt", 0,     0.,   0.,   0.],
         ["Res", 0,     0.,   0.,   ""],
         ["Tot", 0,     0.,   "",   ""]]
# dof row
anova[1][1] = I-1
anova[3][1] = I*J-1
anova[2][1] = anova[3][1] - anova[1][1]
# SS row
anova[2][2] = (J-1)*sum(Syi)
ybar = sum(ybari)/I
for yb in ybari:
    anova[1][2] += (yb - ybar)**2
anova[1][2] *= J
anova[3][2] = anova[1][2] + anova[2][2]
# MS row
anova[1][3] = anova[1][2] / anova[1][1]
anova[2][3] = anova[2][2] / anova[2][1]
# F row
anova[1][4] = anova[1][3] / anova[2][3]

# print ANOVA table col by col to avoid spacing issues
# the damn calculator font isn't unispace anyways
print("\n=== ANOVA ===")
for j in range(1,5):
    for i in range(4):
        print(anova[i][0], anova[i][j])
    print("")

# compute values on the contrasts
a = [[],[],[]] # lin, qdr, cub

print("List the", I, "factors in the linear contrast: ")
while(len(a[0]) < I):
    a[0].append(float(input("")))
print(a[0])
ssa = 0
for ai in a[0]:
    ssa += ai**2
a[0].append(ssa)

print("List the", I, "factors in the quadratic contrast: ")
while(len(a[1]) < I):
    a[1].append(float(input("")))
print(a[1])
ssa = 0
for ai in a[1]:
    ssa += ai**2
a[1].append(ssa)

print("List the", I, "factors in the cubic contrast: ")
while(len(a[2]) < I):
    a[2].append(float(input("")))
print(a[2])
ssa = 0
for ai in a[2]:
    ssa += ai**2
a[2].append(ssa)

# make the contrast table
con = [["   ", "lam", "CI", "SS", "F"],
       ["Lin", 0., "+/- sqrt("+str(I-1)+" * F(alpha,"+str(I-1)+","+str(J*I-4)+")) * ", 0., 0.],
       ["Qua", 0., "+/- sqrt("+str(I-1)+" * F(alpha,"+str(I-1)+","+str(J*I-4)+")) * ", 0., 0.],
       ["Cub", 0., "+/- sqrt("+str(I-1)+" * F(alpha,"+str(I-1)+","+str(J*I-4)+")) * ", 0., 0.]]
# lambdas
for al, aq, ac, yb in zip(a[0], a[1], a[2], ybari):
    con[1][1] += al*yb
    con[2][1] += aq*yb
    con[3][1] += ac*yb
# sum of squares
con[1][3] = J * con[1][1]**2 / a[0][-1]
con[2][3] = J * con[2][1]**2 / a[1][-1]
con[3][3] = J * con[3][1]**2 / a[2][-1]
# confidence intervals
con[1][2] += str(math.sqrt(a[0][-1] * anova[2][3]/ (J)))
con[2][2] += str(math.sqrt(a[1][-1] * anova[2][3]/ (J)))
con[3][2] += str(math.sqrt(a[2][-1] * anova[2][3]/ (J)))
# F
con[1][4] = con[1][3] / anova[2][3]
con[2][4] = con[2][3] / anova[2][3]
con[3][4] = con[3][3] / anova[2][3]

# print col by col
print("\n=== CONTRASTS ===")
for i in range(1, 5):
    for j in range(4):
        print(con[j][0], con[j][i])
    print("")

# compare every mean
print("=== COMPARISONS ===")
diff = [["Comb", "Diff"]]
for i in range(I):
    for j in range(i+1, I):
        this = []
        this.append(str(i+1)+":"+str(j+1))
        this.append(ybari[i]-ybari[j])
        this.append("+/- t(alpha/2,"+str(I*J-1)+") * " + str(math.sqrt(anova[2][3]*2/J)))
        this.append("+/- sqrt("+str(I-1)+" * F(alpha,"+str(I-1)+","+str(J*I-4)+")) * " + str(math.sqrt(2 * anova[2][3]/ (J))))
        diff.append(this)
# do the printing thing
for i in range(1, 2):
    for j in range(len(diff)):
        print(diff[j][0], diff[j][i])
print("LSD: +/- t(alpha/2,"+str(I*J-1)+") *", math.sqrt(anova[2][3]*2/J))
print("Scheffe: +/- sqrt("+str(I-1)+" * F(alpha,"+str(I-1)+","+str(J*I-4)+")) *", math.sqrt(2 * anova[2][3]/ (J)))
print("TW: +/- q("+str(I)+","+str(I*J-I)+") *", math.sqrt(anova[2][3]/J))