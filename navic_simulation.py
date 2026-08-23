import numpy as np
import matplotlib.pyplot as plt

#*--------------------------------------------------------
# GNSS + NavIC POSITION ACCURACY MONTE CARLO SIMULATION
#*--------------------------------------------------------
np.random.seed(42)
NUM_TRIALS = 2000
SPEED_OF_LIGHT = 299_792_458        # m per second

#*-- generate the pseudo range noise standard deviation
#*-- https://insidegnss.com/wp-content/uploads/2018/01/novdec14-UPADHYAY.pdf
def gen_pseudorange(log_cn0):
    cn0 = np.power(10, log_cn0 / 10) # linear signal
    bl = 1.0                        # loop noise bandwidth
    d = 0.5                          # correlator spacing (1 chip)
    t_int = 20.0 / 1000.0            # coherent integration time
    p1 = (bl * d) / cn0
    p2 = 1 + 2 / ( (2 - d) * t_int * cn0)
    return 293 * np.sqrt(p1 * p2)

# Pseudorange measurement noise standard deviation
GPS_SIGNAL = 31
GALILEO_SIGNAL = 28
GLONASS_SIGNAL = 29
BEIDOU_SIGNAL = 30
NAVIC_SIGNAL = 35
SIGMA_GPS = gen_pseudorange(GPS_SIGNAL)       # metres
SIGMA_GALILEO = gen_pseudorange(GALILEO_SIGNAL)   # metres
SIGMA_GLONASS = gen_pseudorange(GLONASS_SIGNAL)   # metres
SIGMA_BEIDOU = gen_pseudorange(BEIDOU_SIGNAL)    # metres
SIGMA_NAVIC = gen_pseudorange(NAVIC_SIGNAL)     # metres

# ------------------------------------------------------------
# Measurement standard deviation
# ------------------------------------------------------------
def get_measurement_sigma(satellite_name):
    if satellite_name.startswith("GPS"):
        return SIGMA_GPS
    elif satellite_name.startswith("GAL"):
        return SIGMA_GALILEO
    elif satellite_name.startswith("GLO"):
        return SIGMA_GLONASS
    elif satellite_name.startswith("BEI"):
            return SIGMA_BEIDOU
    elif satellite_name.startswith("NAVIC"):
        return SIGMA_NAVIC
    else:
        return 3.0

# Number of satellites in each constellation
NUM_GPS = 6
NUM_GALILEO = 4
NUM_GLONASS = 3
NUM_BEIDOU = 9
NUM_NAVIC = 7

#*-------------------------------------------------------------------
# Convert elevation and azimuth into East (x), North (y), and Up (z)
#  ------------------------------------------------------------------
def elevation_azimuth_to_los(elevation_deg, azimuth_deg): 
    el = np.radians(elevation_deg)
    az = np.radians(azimuth_deg)
    east = np.cos(el) * np.sin(az)
    north = np.cos(el) * np.cos(az)
    up = np.sin(el)
    return np.array([east, north, up])

# ------------------------------------------------------------
# Generate a set of GNSS satellites with random elevation and azimuth
# ------------------------------------------------------------
def generate_satellites(num_satellites, name="GNSS"):
    MIN_ELEVATION = 10.0
    MAX_ELEVATION = 89.0
    MIN_AZIMUTH = 0.0
    MAX_AZIMUTH = 360.0
    if name == "NAVIC":     # limit the azimuth to southern sky and elevation
        MIN_ELEVATION = 25.0
        MAX_ELEVATION = 80.0
        MIN_AZIMUTH = 110.0 
        MAX_AZIMUTH = 250.0

    satellites = []
    for i in range(num_satellites):     # all values between min and max have the same prob.
        elevation = np.random.uniform(MIN_ELEVATION, MAX_ELEVATION)
        azimuth = np.random.uniform(MIN_AZIMUTH, MAX_AZIMUTH)
        los = elevation_azimuth_to_los(elevation, azimuth)
        satellites.append({"name": f"{name}_{i + 1}", "elevation": elevation,
            "azimuth": azimuth, "los": los})
    return satellites

# ------------------------------------------------------------
# Create satellite constellations
# ------------------------------------------------------------
gps_satellites = generate_satellites(NUM_GPS, "GPS")
galileo_satellites = generate_satellites(NUM_GALILEO,"GAL")
glonass_satellites = generate_satellites(NUM_GLONASS,"GLO")
beidou_satellites = generate_satellites(NUM_BEIDOU,"BEI")
navic_satellites = generate_satellites(NUM_NAVIC, "NAVIC")

#*------------------------------------------------------------
# Build the GNSS geometry matrix
#*------------------------------------------------------------
def build_geometry_matrix(satellites):
    """
    Construct the geometry matrix G.
    State vector: [East position error, North position error, Up position error, receiver clock bias]
    For each satellite:
        delta_rho = -u_e * delta_E -u_n * delta_N -u_u * delta_U + c * delta_t
    Clock bias is 1 meter
    """
    G = []
    for satellite in satellites:
        u = satellite["los"]
        G.append([-u[0], -u[1], -u[2], 1.0])
    return np.array(G)

# ------------------------------------------------------------
# Weighted DOP calculation
# ------------------------------------------------------------
def calculate_dops(satellites):
    G = build_geometry_matrix(satellites)
    sigmas = np.array([get_measurement_sigma(s["name"]) for s in satellites])
    W = np.diag(1.0 / sigmas**2)    # weight matrix
    try:
        Q = np.linalg.inv(G.T @ W @ G)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan

    # Position covariance terms
    q_ee = Q[0, 0]
    q_nn = Q[1, 1]
    q_uu = Q[2, 2]

    hdop = np.sqrt(q_ee + q_nn)
    vdop = np.sqrt(q_uu)
    pdop = np.sqrt(q_ee + q_nn + q_uu)

    return hdop, vdop, pdop

# ------------------------------------------------------------
# Monte Carlo position solution
# ------------------------------------------------------------

def monte_carlo_accuracy(satellites, num_trials=1000):
    """
    Simulate pseudorange errors and estimate resulting
    receiver position errors using weighted least squares.

    Returns: horizontal RMS error, vertical RMS error, 3D RMS error 95% horizontal error
    """

    G = build_geometry_matrix(satellites)
    sigmas = np.array([get_measurement_sigma(s["name"]) for s in satellites])
    W = np.diag(1.0 / sigmas**2)

    try:
        normal_matrix_inv = np.linalg.inv(G.T @ W @ G)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan, np.nan

    # Weighted least-squares estimator
    H = normal_matrix_inv @ G.T @ W

    horizontal_errors = []
    vertical_errors = []
    error_3d = []

    for _ in range(num_trials):

        # Generate pseudorange measurement noise
        measurement_noise = np.random.normal(loc=0.0, scale=sigmas)

        # Estimated state error: [East, North, Up, Clock]
        state_error = H @ measurement_noise

        east_error = state_error[0]
        north_error = state_error[1]
        up_error = state_error[2]

        horizontal_error = np.sqrt(east_error**2 + north_error**2)
        vertical_error = abs(up_error)
        position_error_3d = np.sqrt(east_error**2 + north_error**2 + up_error**2)

        horizontal_errors.append(horizontal_error)
        vertical_errors.append(vertical_error)
        error_3d.append(position_error_3d)

    horizontal_errors = np.array(horizontal_errors)
    vertical_errors = np.array(vertical_errors)
    error_3d = np.array(error_3d)

    horizontal_rms = np.sqrt(np.mean(horizontal_errors**2))
    vertical_rms = np.sqrt(np.mean(vertical_errors**2))
    error_3d_rms = np.sqrt(np.mean(error_3d**2))
    return (horizontal_rms, vertical_rms, error_3d_rms)

# ------------------------------------------------------------
# Define the scenarios
# ------------------------------------------------------------

other_gnss = (gps_satellites + galileo_satellites + glonass_satellites + beidou_satellites)

scenarios = {
    "Other GNSS only": other_gnss,
    "Other GNSS + 1 NavIC": other_gnss + navic_satellites[:1],
    "Other GNSS + 2 NavIC": other_gnss + navic_satellites[:2],
    "Other GNSS + 3 NavIC": other_gnss + navic_satellites[:3],
    "Other GNSS + 4 NavIC": other_gnss + navic_satellites[:4],
    "Other GNSS + 5 NavIC": other_gnss + navic_satellites[:5],
    "Other GNSS + 6 NavIC": other_gnss + navic_satellites[:6],
    "Other GNSS + 7 NavIC": other_gnss + navic_satellites[:7]
}

# ------------------------------------------------------------
# Run the simulation
# ------------------------------------------------------------
results = {}
print()
print("=" * 72)
print("GNSS + NAVIC POSITION ACCURACY SIMULATION")
print("=" * 72)

for scenario_name, satellites in scenarios.items():
    hdop, vdop, pdop = calculate_dops(satellites)
    (horizontal_rms, vertical_rms, error_3d_rms) = \
        monte_carlo_accuracy(satellites, NUM_TRIALS)

    results[scenario_name] = {
        "num_satellites": len(satellites), "hdop": hdop, "vdop": vdop, "pdop": pdop,
        "horizontal_rms": horizontal_rms, "vertical_rms": vertical_rms, "error_3d_rms": error_3d_rms}

    print(f"\n{scenario_name}")
    print("-" * 72)
    print(f"Number of satellites : {len(satellites)}")
    print(f"HDOP                 : {hdop:.3f}")
    print(f"VDOP                 : {vdop:.3f}")
    print(f"PDOP                 : {pdop:.3f}")
    print(f"Horizontal RMS error : {horizontal_rms:.3f} m")
    print(f"Vertical RMS error   : {vertical_rms:.3f} m")
    print(f"3D RMS error         : {error_3d_rms:.3f} m")

# ------------------------------------------------------------
# Calculate percentage improvement
# ------------------------------------------------------------

baseline = results["Other GNSS only"]
print()
print("=" * 72)
print("IMPROVEMENT RELATIVE TO OTHER GNSS ONLY")
print("=" * 72)
DIR = "/home/mkonchady/books/gps/ch10/fig6/"
f = open(DIR + "plot6.dat", "w")
sat_num = 1
for scenario_name in scenarios:
    if scenario_name == "Other GNSS only":
        continue
    result = results[scenario_name]
    horizontal_improvement = (100.0 * (baseline["horizontal_rms"] - result["horizontal_rms"])
        / baseline["horizontal_rms"])

    vertical_improvement = (100.0 * (baseline["vertical_rms"] - result["vertical_rms"])
        / baseline["vertical_rms"])

    error_3d_improvement = (100.0 * (baseline["error_3d_rms"] - result["error_3d_rms"])
        / baseline["error_3d_rms"])

    print(f"\n{scenario_name}")
    print(f"Horizontal RMS improvement: " f"{horizontal_improvement:.1f}%")
    print(f"Vertical RMS improvement  : " f"{vertical_improvement:.1f}%")
    print(f"3D RMS improvement        : " f"{error_3d_improvement:.1f}%")
    f.write(str(sat_num) + "," + str(horizontal_improvement) + "," +  str(error_3d_improvement) + "\n")
    sat_num = sat_num + 1
f.close()
print ("done")