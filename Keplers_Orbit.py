import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, pi, sqrt
from scipy.integrate import quad

# fmt: off
u = 0.0  # x-position of the center
v = 0  # y-position of the center
a = 10.0  # radius on the x-axis
b = 1  # radius on the y-axis

t = np.linspace(0, 2 * pi, 100)
plt.plot(u + a * np.cos(t), v + b * np.sin(t))
plt.grid(color="lightgray", linestyle="--")
plt.plot(-5, 0, marker="o", markersize=50, markeredgecolor="orange", markerfacecolor="yellow")
plt.text(-5.6, -0.04, "Sun", fontweight="bold")
plt.plot(-10, 0, marker="o", markersize=30, markeredgecolor="red", markerfacecolor="red")
plt.text(-11, -0.25, "perihelion", fontweight="bold", fontsize=12)
plt.plot(10, 0, marker="o", markersize=30, markeredgecolor="green", markerfacecolor="green")
plt.text(6.5, -0.25, "aphelion", fontweight="bold", fontsize=12)
plt.plot(-5, 0.875,
    marker="o",
    markersize=10,
    markeredgecolor="black",
    markerfacecolor="grey",
)
plt.text(-5, 0.60, "Average day, between \n perihelion/aphelion")
point1 = [-6.6, -0]
point2 = [-9.9, -0.17]
x_values = [point1[0], point2[0]]
y_values = [point1[1], point2[1]]
point4 = [-9.9, 0.17]
x_values2 = [point1[0], point4[0]]
y_values2 = [point1[1], point4[1]]
point5 = [-3.4, -0]
point6 = [9.9, -0.17]
x_values3 = [point5[0], point6[0]]
y_values3 = [point5[1], point6[1]]
point7 = [9.9, 0.16]
x_values4 = [point5[0], point7[0]]
y_values4 = [point5[1], point7[1]]
plt.plot(x_values, y_values, color="purple", linewidth=3)
plt.plot(x_values2, y_values2, color="purple", linewidth=3)
plt.plot(x_values3, y_values3, color="orange", linewidth=3)
plt.plot(x_values4, y_values4, color="orange", linewidth=3)
plt.axis("off")
plt.show()

# according to NASA--> https://science.nasa.gov/science-news/science-at-nasa/2001/ast04jan_1 and https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html:
class EarthParameters:
    AU = 149597870600.91  # length of 1 astronomical unit in meters
    peri = 0.9832899 * AU  # perihelion in meters, when earth is closest to the sun
    aph = 1.0167103 * AU  # aphelion in meters, when earth is farthest from the sun
    e = 0.016710219  # eccentricity of earth's orbit
    peri_v = 30288  # perihelion velocity in m/s
    aph_v = 29292  # aphelion veliocity in m/s
    avg_v = 29780
    angular_peri_v = peri_v / peri  # angular velocity at perihelion PER SECOND
    angular_aph_v = aph_v / aph  # angular velocity at aphelion PER SECOND
    angular_avg_v = avg_v / AU  # average angular velocity PER SECOND
    semi_maj = 1.0000001124 * AU  # semimajor axis
    semi_min = 0.9998604869 * AU  # semiminor axis
    orbit = 365.25  # exact number of days to revolve around sun
    T = 365.25 * 86400  # seconds per year
    earth_rotation = 7.292115e-5  # earth's rotational speed in radians/second, according to Wikipedia
    deg_peri = (
        angular_peri_v
    ) * 86400  # rotates through this many radians at perihelion day
    deg_aph = (
        angular_aph_v
    ) * 86400  # rotates through this many radians at aphelion day
    deg_avg = (angular_avg_v) * 86400  # rotates through this many radians on avg day

# according to the conservation of angular momentum, we know that m*v*r is a constant. The mass of the Earth is not changing, so we can say that v*r is a constant. To find the approximate value of this, we use the average orbital speed of the earth and the average distance, which is 1 AU.
constant_Vr = EarthParameters.avg_v * EarthParameters.AU

distance_p_to_a = (EarthParameters.aph - EarthParameters.peri) / 183
# This is the approximate change in distance per day between THIS YEAR’S perihelion (January 4th) and aphelion (July 6th), which is 183 days.

distance_a_to_p = (EarthParameters.aph - EarthParameters.peri) / 180
# This is the approximate change in distance per day between THIS YEAR’S aphelion (July 6th) and NEXT YEAR’S perihelion (January 2nd), which is 180 days.

previous_aph_to_current_peri_distance_per_day = (EarthParameters.aph - EarthParameters.peri) / 184
# This is the approximate change in distance per day between LAST YEAR’S aphelion (July 4th) and THIS YEAR’S perihelion (January 4th), which is 184 days.

previous_aph_to_current_peri_radians_per_day = (pi / 184)  # The earth will cover pi radians around the ellipse from aphelion of last year to perihelion of this year; this metric gives the average number of radians it travels per day.

df = pd.DataFrame(columns=["Date", "Earth_Orbital_Velocity", "Distance"])

df["Date"] = pd.date_range(start="1/1/2023", end="12/31/2023", freq="D")  # creating df column that includes all 365 days of the year
df["Date"] = pd.to_datetime(df["Date"])
df["Date"] = df["Date"].dt.strftime("%Y-%m-%d")

df["Distance"][3] = EarthParameters.peri  # filling in perihelion distance
df["Distance"][186] = EarthParameters.aph  # filling in aphelion distance
# #here, we fill in the approximate distance from the Sun in meters for the first three days of the year, before perihelion:
df["Distance"][2] = EarthParameters.peri + previous_aph_to_current_peri_distance_per_day
df["Distance"][1] = (EarthParameters.peri + 2 * previous_aph_to_current_peri_distance_per_day)
df["Distance"][0] = (EarthParameters.peri + 3 * previous_aph_to_current_peri_distance_per_day)

# here, we fill in approximate values of how the distances changes per day between the Earth and the Sun.
for i in range(4, 186):
    df["Distance"][i] = df["Distance"][i - 1] + distance_p_to_a
for i in range(187, 365):
    df["Distance"][i] = df["Distance"][i - 1] - distance_a_to_p

df["Earth_Orbital_Velocity"] = (constant_Vr / df["Distance"])  # this is earth's orbital velocity. Note that this is different than rotation velocity (spinning on its own axis)
df["Angular_Velocity"] = (df["Earth_Orbital_Velocity"] / df["Distance"])  # this is the angular velocity in radians/second
df["Radians_of_Rotation"] = (df["Angular_Velocity"] * 86400)  # radians the Earth rotates through each day based on our definition of a day being 86400 seconds.

# https://www.johndcook.com/blog/2022/11/02/elliptic-arc-length/
def arc_length_in_meters(T0, T1, a, b):
    f = lambda t: sqrt((a**2 * sin(t) ** 2 + b**2 * cos(t) ** 2))
    return quad(f, T0, T1)

arc_lengths = []
for i in df.iloc:
    start = 0  # we can start at 0 because the only difference in angle we care about is the radians of rotation. This makes the integral (and the code) simpler!
    end = i.Radians_of_Rotation
    X = arc_length_in_meters(start, end, EarthParameters.semi_maj, EarthParameters.semi_min)
    X = X[0]
    arc_lengths.append(X)
df["Arc_Length"] = arc_lengths

df["Rotation_Time"] = df["Arc_Length"] / df["Earth_Orbital_Velocity"]
df["Difference_from_defined_Day_in_Minutes"] = (df["Rotation_Time"] - 86400) / 60
df["Difference_per_hour_in_minutes"] = df["Difference_from_defined_Day_in_Minutes"] / 24
df["Difference_per_minute"] = df["Difference_per_hour_in_minutes"] / 60
df["Difference_per_second"] = df["Difference_per_minute"] / 60


