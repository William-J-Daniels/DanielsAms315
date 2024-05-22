import math

# get info from problem
I = int(input("Number of levels: "))
print(I)

J = int(input("Observations per level: "))
print(J)

ybarn = []
Syn = []
for i in range(I):
    ybarn.append(float(input("ybar"+str(i+1)+": ")))
    print(ybarn[-1])

for i in range(I):
    Syn.append(float(input("Sy"+str(i+1)+": ")))
    print(Syn[-1])

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
anova[2][2] = (J-1) * sum(Syn)
ybar = sum(ybarn) / I
for yb in ybarn:
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
print("=== ANOVA ===")
for j in range(1,5):
    for i in range(4):
        print(anova[i][0], anova[i][j])
    print("")


# fishers least significant difference
print("=== FISHERS LSD ===")
fld = [["Comb", "Diff"]]
for i in range(I):
    for j in range(i+1, I):
        fld.append([str(i+1)+":"+str(j+1), (ybarn[i]-ybarn[j])])
for row in fld:
    print(row)
print("Diff +/- t(alpha/2,"+str(I*J-1)+") *", math.sqrt(anova[2][3]*2/J))