import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import fnmatch
import pygmt
from matplotlib.colors import ListedColormap
from scipy.ndimage import label
import pandas as pd
from scipy.stats import linregress
import math
from matplotlib.colors import to_rgba
import copy
from scipy.interpolate import interp1d
from matplotlib.backends.backend_pdf import PdfPages
import glob
from scipy.stats import pearsonr



def seconds_to_years(year_n):
    seconds_per_year = year_n / (365.25 * 24 * 3600)
    return seconds_per_year

def find_events(masked_grid, DX, N2, plotYN):
    binary_mask = ~np.isnan(masked_grid)
    labeled_fault_matrix, ngroups = label(binary_mask)
    timestep_event = []
    rupture_length_pixels = []  # length of the rupture in pixels (xaxis)
    group_id = []

    for groupi in range(1, ngroups + 1):
        group_indices = np.argwhere(labeled_fault_matrix == groupi)
        timestep = group_indices[0, 0]  # onset of rupture
        x_indices = group_indices[:, 1]
        rupture_length = (
            np.ptp(x_indices) + 1
        )  # ptp calculates the range, +1 to include both endpoints
        timestep_event.append(timestep)
        rupture_length_pixels.append(rupture_length)
        group_id.append(groupi)

    total_npixels = masked_grid.shape[1]
    pixel_size = (DX * N2) / total_npixels  # size of pixels

    if plotYN == "Yes":
        plt.figure(dpi=300)
        cmap = plt.get_cmap("tab20", ngroups + 1)  # Add 1 to accommodate event 0
        cmap.colors[0] = to_rgba("white")  #
        plt.imshow(
            np.flipud(labeled_fault_matrix),
            cmap=cmap,
            aspect="auto",
            interpolation="none",
        )
        plt.colorbar(label="Group ID")
        plt.xlabel("Position")
        plt.ylabel("Time step #")
        plt.xticks(rotation=45)
        plt.gca().set_facecolor("white")

    return (timestep_event, rupture_length_pixels, group_id, total_npixels, pixel_size)


def potency(grid, DX):
    # grid = np.flipud(grid)
    potency = []
    for row in grid:
        dsum = np.sum(row * DX)
        potency.append(dsum)
    return potency


def estimate_Ru(shear_mod, Dc, b, a, sigma_o, W):
    h_star = (shear_mod * Dc) / ((b - a) * sigma_o)
    Ru = W / h_star
    return Ru


def estimate_h_star(shear_mod, Dc, b, a, sigma_o):
    h_star = (shear_mod * Dc) / ((b - a) * sigma_o)
    return h_star


def estimate_D_W_ratio(D, W):
    D_W = D / W
    return D_W


def plot_interevent_vs_rupture(fault_data):
    grouped = fault_data.groupby(["D", "Dc"])
    results = []
    plt.figure(dpi=300)
    for (D_value, Dc_value), group in grouped:
        group = group[group["Event_rupture_length_km"] >= 0.5]
        group = group.sort_values(by="Event_time_seconds")
        group["Interevent_time"] = (
            group["Event_time_seconds"].diff().shift(-1)
        )  #  time to the next rupture
        # group = group[group['Interevent_time'] >= 10000]
        marker = marker_shapes.get(D_value, "o")
        color = dc_colors.get(Dc_value, "black")
        alpha = 0.1 if Dc_value in [0.001, 0.003] else 0.5
        plt.scatter(
            group["Event_rupture_length_km"],
            group["Interevent_time"] / 365 / 24 / 60 / 60,
            alpha=alpha,
            color=color,
            marker=marker,
            label=f"D={D_value}, Dc={Dc_value}",
        )
        group_for_export = group[
            ["Event_rupture_length_km", "Interevent_time"]
        ].dropna()
        group_for_export["Dc"] = Dc_value
        group_for_export["D"] = D_value
        results.append(group_for_export)

    plt.ylabel("Time to next rupture (years)")
    plt.xlabel("Rupture length (km)")

    results_df = pd.concat(results).reset_index(drop=True)
    return results_df


def fit_and_plot_combined_line(fault1_times_length, fault2_times_length):
    fault1_times_length["Interevent_time_years"] = (
        fault1_times_length["Interevent_time"] / 365 / 24 / 60 / 60
    )
    fault2_times_length["Interevent_time_years"] = (
        fault2_times_length["Interevent_time"] / 365 / 24 / 60 / 60
    )

    combined_data = pd.concat(
        [fault1_times_length, fault2_times_length], ignore_index=True
    )

    plt.figure(dpi=300)

    grouped = combined_data.groupby(["D", "Dc"])

    for (D_value, Dc_value), group in grouped:
        color = dc_colors.get(Dc_value, "black")
        marker = marker_shapes.get(D_value, "o")
        alpha = 0.3 if Dc_value in [0.001, 0.003] else 0.3

        plt.scatter(
            group["Event_rupture_length_km"],
            group["Interevent_time_years"],
            alpha=alpha,
            color=color,
            marker=marker,
            edgecolor="none",
            label=f"D={D_value}, Dc={Dc_value}",
        )

    slope, intercept, r_value, p_value, std_err = linregress(
        combined_data["Event_rupture_length_km"], combined_data["Interevent_time_years"]
    )

    xfit = np.linspace(
        min(combined_data["Event_rupture_length_km"]),
        max(combined_data["Event_rupture_length_km"]),
        1000,
    )
    fit_line = slope * xfit + intercept
    plt.plot(
        xfit,
        fit_line,
        label=f"Fit: y={slope:.2f}x + {intercept:.2f}",
        color="slategray",
    )
    plt.xlabel("Rupture length (km)")
    plt.ylabel("Interevent time (years)")

    shape_legend = [
        plt.Line2D(
            [0], [0], color="black", marker=marker_shapes[d], linestyle="", markersize=6
        )
        for d in marker_shapes.keys()
    ]
    shape_legend_labels = [f"{d} km" for d in marker_shapes.keys()]
    color_legend = [
        plt.Line2D(
            [0], [0], color=dc_colors[dc], marker="o", linestyle="", markersize=6
        )
        for dc in dc_colors.keys()
    ]
    color_legend_labels = [f"{dc} m" for dc in dc_colors.keys()]

    legend1 = plt.legend(
        shape_legend,
        shape_legend_labels,
        fontsize=6,
        title="D",
        loc="upper right",
        ncol=2,
    )
    plt.gca().add_artist(legend1)

    legend2 = plt.legend(
        color_legend, color_legend_labels, fontsize=6, title=r"$D_c$", loc="upper left"
    )
    plt.gca().add_artist(legend2)
    return slope, intercept, r_value


def calculate_h_2D(G, L, b, a, sigma):
    h_2D_star = (2 * G * L * b) / (
        math.pi * sigma * (b - a) ** 2
    )  # Rubin and Ampuero 2D 2005
    return h_2D_star


def cycle_duration(b, a, sigma, G, Vpl, Lweak, L, D, Vdyn):
    t = (
        (((b - a) * sigma) / (G * Vpl))
        * np.sqrt((Lweak - (Lweak - calculate_h_2D(G, L, b, a, sigma))) * D)
        * np.log(Vdyn / Vpl)
    )
    return t


def cycles_colormap(m=None):
    """
    Original code by Sylvain Barbot, ceded by Baoning Wu.
    Generates a colormap for representing velocity during seismic cycles.
    """
    if m is None:
        m = plt.rcParams["image.lut"]
    cpt = np.array(
        [
            [-12.0, 0, 0, 0],
            [-11.0, 0, 0, 0],
            [-9.3, 106, 135, 196],
            [-8.8, 135, 164, 224],
            [-3.0, 247, 236, 44],
            [-1.0, 239, 64, 35],
            [0.0, 128, 21, 23],
            [1.0, 50, 21, 23],
        ]
    )

    x = -12 + np.linspace(0, 1, m) * 13.0

    # interp RGB channels
    r_interp = interp1d(cpt[:, 0], cpt[:, 1] / 255, kind="linear")
    g_interp = interp1d(cpt[:, 0], cpt[:, 2] / 255, kind="linear")
    b_interp = interp1d(cpt[:, 0], cpt[:, 3] / 255, kind="linear")

    r = r_interp(x)
    g = g_interp(x)
    b = b_interp(x)
    colormap = np.stack([r, g, b], axis=1)

    return colormap


### KDC phase code
def get_phase_time_series(tq, t):
    return np.array([get_phase(itq, t) for itq in tq])


def catalog2catalog_order_parameter(cq, c):
    phases = get_phase_time_series(cq.catalog.time.values, c.catalog.time.values)
    phases = phases[~np.isnan(phases)]
    if np.any(phases):
        order_parameter = get_order_parameter(phases)
    else:
        order_parameter = np.nan

    return order_parameter


def catalog2catalog_phase(cq, c, remove_nan=True):
    phases = get_phase_time_series(cq.catalog.time.values, c.catalog.time.values)
    if remove_nan:
        phases = phases[~np.isnan(phases)]
    return phases


def phase_analysis(repeaters, earthquakes, search_radius=0.1, inner_radius=0):

    order_parameter = []
    number_of_neighboring_repeaters = []
    number_of_neighboring_earthquakes = []
    earthquake_order_parameter = []
    raw_phases = []
    raw_times = []
    phases = []
    earthquake_phases = []
    delta_omega = []

    families = repeaters.get_families()

    for family in families:
        family_ID = family.catalog.family.values[0]

        earthquakes.catalog["distance"] = family.get_nearest_neighbor_distance(
            earthquakes.catalog
        )

        repeaters.catalog["distance"] = family.get_nearest_neighbor_distance(
            repeaters.catalog
        )

        # make a deep copy of the catalog
        neighboring_earthquakes = copy.deepcopy(earthquakes)
        neighboring_repeaters = copy.deepcopy(repeaters)

        neighboring_repeaters.catalog = neighboring_repeaters.catalog.loc[
            (neighboring_repeaters.catalog.distance < inner_radius + search_radius)
            & (neighboring_repeaters.catalog.distance > inner_radius)
        ]

        neighboring_earthquakes.catalog = neighboring_earthquakes.catalog.loc[
            (neighboring_earthquakes.catalog.distance < inner_radius + search_radius)
            & (neighboring_earthquakes.catalog.distance > inner_radius)
        ]

        neighboring_repeaters.catalog = neighboring_repeaters.catalog.loc[
            neighboring_repeaters.catalog.family != family_ID
        ]

        order_parameter.append(
            catalog2catalog_order_parameter(neighboring_repeaters, family)
        )

        earthquake_order_parameter.append(
            catalog2catalog_order_parameter(neighboring_earthquakes, family)
        )

        number_of_neighboring_repeaters.append(len(neighboring_repeaters))
        number_of_neighboring_earthquakes.append(len(neighboring_earthquakes))

        phases.append(catalog2catalog_phase(neighboring_repeaters, family))
        raw_phases.append(
            catalog2catalog_phase(neighboring_repeaters, family, remove_nan=False)
        )
        raw_times.append(neighboring_repeaters.catalog.time.values)
        earthquake_phases.append(catalog2catalog_phase(neighboring_earthquakes, family))

        delta_omega.append(
            (
                1 / family.catalog.RCm.values[0]
                - 1 / neighboring_repeaters.catalog.RCm.values
            )
            * family.catalog.RCm.values[0]
        )

    order_parameter = np.array(order_parameter)
    earthquake_order_parameter = np.array(earthquake_order_parameter)
    raw_phases = np.concatenate(raw_phases)
    raw_times = np.concatenate(raw_times)
    earthquake_phases = np.concatenate(earthquake_phases)

    return (
        order_parameter,
        earthquake_order_parameter,
        phases,
        raw_phases,
        raw_times,
        earthquake_phases,
        number_of_neighboring_repeaters,
        number_of_neighboring_earthquakes,
        delta_omega,
    )


def get_phase(tq, t: np.ndarray):
    # note: not vectorized

    if len(t) < 2:
        return np.nan

    if tq < min(t) or tq > max(t):
        return np.nan

    dt = t - tq

    if np.min(np.abs(dt)) == 0:
        index = np.argmin(np.abs(dt))
        t1 = 0
        if index == 0:
            T = dt[index + 1]
        elif index == len(dt) - 1:
            T = dt[index - 1]
        else:
            T = np.mean(dt[[index - 1, index + 1]])

    else:
        dt_pos = dt[dt >= 0]
        dt_neg = dt[dt < 0]

        t1 = min(-dt_neg)
        t2 = min(dt_pos)

        T = t1 + t2

    phase = 2 * np.pi * t1 / T

    return phase


def get_order_parameter(phases):
    # Returns alignment coefficient = mean(cos(θ)), NOT the Kuramoto order parameter |mean(exp(iθ))|
    order_real = np.nanmean(np.cos(phases))
    #  np.abs(np.nanmean(np.exp(1j * phases)))
    return order_real


def measure_intevent_time_f(times):
    times = times
    sorted_times = np.sort(times)
    interevent_times = np.diff(sorted_times)
    return interevent_times


def get_phase_unbounded(tq, t: np.ndarray):
    dt = t - tq
    closest_idx = np.argmin(np.abs(dt))  # closest event in time sequence
    cycle_period = np.mean(np.diff(t))
    # calc how many cycles away from the closest event
    cycles_away = dt[closest_idx] / cycle_period
    phase = (cycles_away % 1) * 2 * np.pi

    return phase


def markdown_table_to_pdf(markdown_table, output_pdf="table.pdf"):
    # markdown table to pandas df
    # First, split the table into lines and remove empty lines
    lines = [
        line.strip() for line in markdown_table.strip().split("\n") if line.strip()
    ]
    lines = [line for line in lines if not all(c in "-|" for c in line)]

    data = []
    for line in lines:
        cells = [cell.strip() for cell in line.split("|") if cell.strip()]
        data.append(cells)

    df = pd.DataFrame(data[1:], columns=data[0])

    fig, ax = plt.subplots(figsize=(10, len(df) * 0.5 + 1))

    ax.axis("tight")
    ax.axis("off")

    table = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        cellLoc="center",
        loc="center",
        bbox=[0, 0, 1, 1],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 1.5)
    with PdfPages(output_pdf) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close()


def Griffin_plot_interevent_times(interevent_times_df, ax, marker_shapes, dc_colors):
    # for each combination of Dc and D
    grouped = interevent_times_df.groupby(["Dc", "D"])
    for (Dc, D), group in grouped:

        def process_times(times):
            if isinstance(times, str):
                return [float(x) for x in times.split(",")]
            elif isinstance(times, (list, np.ndarray)):
                return [float(x) for x in times]
            else:
                return [float(times)]  # Handle single value case

        times_f1 = process_times(group["Inter-event times fault 1 (years)"].iloc[0])
        times_f2 = process_times(group["Inter-event times fault 2 (years)"].iloc[0])

        N_f1 = len(times_f1)
        mu_1_f1 = np.mean(times_f1[:-1])
        mu_2_f1 = np.mean(times_f1[1:])
        sigma_1_f1 = np.std(times_f1[:-1])
        sigma_2_f1 = np.std(times_f1[1:])
        numerator_f1 = (times_f1[:-1] - mu_1_f1) * (times_f1[1:] - mu_2_f1)
        memory_f1 = np.sum(numerator_f1) / ((N_f1 - 1) * sigma_1_f1 * sigma_2_f1)
        # burstiness
        burstiness_f1 = (np.std(times_f1) - np.mean(times_f1)) / (
            np.std(times_f1) + np.mean(times_f1)
        )

        N_f2 = len(times_f2)
        mu_1_f2 = np.mean(times_f2[:-1])
        mu_2_f2 = np.mean(times_f2[1:])
        sigma_1_f2 = np.std(times_f2[:-1])
        sigma_2_f2 = np.std(times_f2[1:])
        numerator_f2 = (times_f2[:-1] - mu_1_f2) * (times_f2[1:] - mu_2_f2)
        memory_f2 = np.sum(numerator_f2) / ((N_f2 - 1) * sigma_1_f2 * sigma_2_f2)
        burstiness_f2 = (np.std(times_f2) - np.mean(times_f2)) / (
            np.std(times_f2) + np.mean(times_f2)
        )

        marker_shape = marker_shapes.get(D, "o")
        color = dc_colors.get(Dc, "black")
        ax.scatter(
            memory_f1,
            burstiness_f1,
            color=color,
            marker=marker_shape,
            s=20,
            alpha=0.6,
            edgecolor="none",
        )
        ax.scatter(
            memory_f2,
            burstiness_f2,
            color=color,
            marker=marker_shape,
            s=20,
            alpha=0.6,
            edgecolor="none",
        )


### functions for processing paleoseismic data


def get_moment_files_for_pair(fault1, fault2, moment_dir, dl_dir):
    result = {
        "fault1_file": None,
        "fault2_file": None,
        "dl_file": None,
        "fault1_name": fault1,
        "fault2_name": fault2,
    }

    for csv_file in glob.glob(os.path.join(moment_dir, "*.csv")):
        basename = os.path.basename(csv_file)
        if fault1 in basename:
            result["fault1_file"] = csv_file
        if fault2 in basename:
            result["fault2_file"] = csv_file

    for xlsx_file in glob.glob(os.path.join(dl_dir, "*.xlsx")):
        basename = os.path.basename(xlsx_file)
        if fault1 in basename and fault2 in basename:
            result["dl_file"] = xlsx_file
            break

    return result


# MC for data corr coeff and alignment uncertainties for paleoseismic data
def alignment_MC(fault1_tmin, fault1_tmax, fault2_tmin, fault2_tmax, n_mc=1000):
    f1_mid = np.sort((fault1_tmin + fault1_tmax) / 2)
    f2_mid = np.sort((fault2_tmin + fault2_tmax) / 2)

    phases_12_mid = np.array([get_phase(t, f1_mid) for t in f2_mid])
    phases_12_mid = phases_12_mid[~np.isnan(phases_12_mid)]
    ac_12_mid = np.nanmean(np.cos(phases_12_mid)) if len(phases_12_mid) > 0 else np.nan

    phases_21_mid = np.array([get_phase(t, f2_mid) for t in f1_mid])
    phases_21_mid = phases_21_mid[~np.isnan(phases_21_mid)]
    ac_21_mid = np.nanmean(np.cos(phases_21_mid)) if len(phases_21_mid) > 0 else np.nan

    ac_12_samples = []
    ac_21_samples = []
    for _ in range(n_mc):
        f1_times = np.sort(np.random.uniform(fault1_tmin, fault1_tmax))
        f2_times = np.sort(np.random.uniform(fault2_tmin, fault2_tmax))
        phases_12 = np.array([get_phase(t, f1_times) for t in f2_times])
        phases_12 = phases_12[~np.isnan(phases_12)]
        phases_21 = np.array([get_phase(t, f2_times) for t in f1_times])
        phases_21 = phases_21[~np.isnan(phases_21)]
        if len(phases_12) > 0:
            ac_12_samples.append(np.nanmean(np.cos(phases_12)))
        if len(phases_21) > 0:
            ac_21_samples.append(np.nanmean(np.cos(phases_21)))

    def _summarise(samples, midpoint):
        if len(samples) == 0:
            return dict(midpoint=np.nan, p5=np.nan, p95=np.nan, n_valid=0)
        return {
            "midpoint": midpoint,
            "p5": np.percentile(samples, 5),
            "p95": np.percentile(samples, 95),
            "n_valid": len(samples),
        }

    return {
        "ac_12": _summarise(ac_12_samples, ac_12_mid),
        "ac_21": _summarise(ac_21_samples, ac_21_mid),
    }


def slip_rate_correlation_MC(f1_upper, f1_lower, f2_upper, f2_lower, n_mc=1000):
    times = np.union1d(
        np.union1d(f1_upper["time"], f1_lower["time"]),
        np.union1d(f2_upper["time"], f2_lower["time"]),
    )

    f1u_sort = np.argsort(f1_upper["time"])
    f1l_sort = np.argsort(f1_lower["time"])
    f2u_sort = np.argsort(f2_upper["time"])
    f2l_sort = np.argsort(f2_lower["time"])

    f1_hi = np.interp(times, f1_upper["time"][f1u_sort], f1_upper["slip_rate"][f1u_sort])
    f1_lo = np.interp(times, f1_lower["time"][f1l_sort], f1_lower["slip_rate"][f1l_sort])
    f2_hi = np.interp(times, f2_upper["time"][f2u_sort], f2_upper["slip_rate"][f2u_sort])
    f2_lo = np.interp(times, f2_lower["time"][f2l_sort], f2_lower["slip_rate"][f2l_sort])

    f1_lo, f1_hi = np.minimum(f1_lo, f1_hi), np.maximum(f1_lo, f1_hi)
    f2_lo, f2_hi = np.minimum(f2_lo, f2_hi), np.maximum(f2_lo, f2_hi)

    corr_samples = []

    for _ in range(n_mc):
        s1 = f1_lo + np.random.uniform(0, 1, len(times)) * (f1_hi - f1_lo)
        s2 = f2_lo + np.random.uniform(0, 1, len(times)) * (f2_hi - f2_lo)

        if np.std(s1) == 0 or np.std(s2) == 0:
            continue

        r, _ = pearsonr(s1, s2)
        corr_samples.append(r)

    if len(corr_samples) == 0:
        return dict(corr_mean=np.nan, corr_std=np.nan, corr_p5=np.nan, corr_p95=np.nan)

    corr_samples = np.array(corr_samples)
    return {
        "corr_mean": np.mean(corr_samples),
        "corr_std":  np.std(corr_samples),
        "corr_p5":   np.percentile(corr_samples, 5),
        "corr_p95":  np.percentile(corr_samples, 95),
    }


def moment_correlation_MC(
    f1_years, f1_rates, f2_years, f2_rates, uncertainty_fraction=0.20, n_mc=1000
):
    f1_sort = np.argsort(f1_years)
    f2_sort = np.argsort(f2_years)

    years_common = np.union1d(f1_years, f2_years)
    f1_interp = np.interp(years_common, f1_years[f1_sort], f1_rates[f1_sort])
    f2_interp = np.interp(years_common, f2_years[f2_sort], f2_rates[f2_sort])

    f1_sigma = np.abs(f1_interp) * uncertainty_fraction
    f2_sigma = np.abs(f2_interp) * uncertainty_fraction

    corr_samples = []

    for _ in range(n_mc):
        s1 = f1_interp + np.random.normal(0, f1_sigma)
        s2 = f2_interp + np.random.normal(0, f2_sigma)
        if np.std(s1) == 0 or np.std(s2) == 0:
            continue
        r, _ = pearsonr(s1, s2)
        corr_samples.append(r)

    if len(corr_samples) == 0:
        return dict(corr_mean=np.nan, corr_std=np.nan, corr_p5=np.nan, corr_p95=np.nan)

    corr_samples = np.array(corr_samples)
    return {
        "corr_mean": np.mean(corr_samples),
        "corr_std":  np.std(corr_samples),
        "corr_p5":   np.percentile(corr_samples, 5),
        "corr_p95":  np.percentile(corr_samples, 95),
    }