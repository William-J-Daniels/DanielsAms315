import math

# get inputs from problem
n = int(input("n: "))
print(n)

ybar = float(input("ybar: "))
print(ybar)

Sy = float(input("Sy: "))
print(Sy)

xbar = float(input("xbar: "))
print(xbar)

Sx = float(input("Sx: "))
print(Sx)

r = float(input("r: "))
print(r)


# compute ANOVA values
anova = [["   ", "dof", "SS", "MS", "F"],
         ["Reg", 1,     0.,   0.,   0.],
         ["Res", n-2,   0.,   0.,   ""],
         ["Tot", n-1,   0.,   "",   ""]]
# SS row
anova[3][2] = anova[3][1] * Sy**2
anova[1][2] = anova[3][2] * r**2
anova[2][2] = anova[3][2] - anova[1][2]
# MS row
anova[1][3] = anova[1][2]
anova[2][3] = anova[2][2] / anova[2][1]
#F row
anova[1][4] = anova[1][3] / anova[2][3]


# print ANOVA table col by col to avoid spacing issues
# the damn calculator font isn't unispace anyways
print("=== ANOVA ===")
for j in range(1,5):
    for i in range(4):
        print(anova[i][0], anova[i][j])
    print("")


# compute regression parameters
B1 = r * Sy / Sx
B0 = ybar - B1*xbar
print("=== REGRESSION ===")
print("yhat(x) =", B0, "+", B1, "* x")
print("B1 CI =", B1, "+/- t(alpha/2,"+str(n-1)+") *", math.sqrt(anova[2][3] / (anova[3][1] * Sx**2)))
print("")

print("=== INTERVALS ===")
# CI and PI at specified x
x = float(input("Value of x to find CI and PI: "))
print(x)

print(B0 + B1*x)
print("CI: +/- t(alpha/2,"+str(n-1)+") *", math.sqrt(anova[2][3] * (1/n + (x - xbar)**2 / (anova[3][1] * Sx**2))))
print("PI: +/- t(alpha/2,"+str(n-1)+") *", math.sqrt(anova[2][3] * (1 + 1/n + (x - xbar)**2 / (anova[3][1] * Sx**2))))