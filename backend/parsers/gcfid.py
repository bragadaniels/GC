import io
import numpy as np
from scipy.io import netcdf_file


def _get_nc_attr(var, attr_name, default=None):
    """Extrai atributos como scale_factor e add_offset com segurança."""
    if hasattr(var, attr_name):
        return getattr(var, attr_name)
    if hasattr(var, "_attributes") and attr_name in var._attributes:
        return var._attributes[attr_name]
    return default


def read_gcfid_cdf(file_bytes: bytes, filename: str) -> dict:
    """
    Lê o arquivo .CDF de GC-FID aplicando o scale_factor do detector NetCDF.
    """
    try:
        nc = netcdf_file(io.BytesIO(file_bytes), "r", mmap=False)
        var_ord = nc.variables["ordinate_values"]

        # Fator de escala do detector (vital para bater com a área do software do aparelho!)
        scale_factor = float(_get_nc_attr(var_ord, "scale_factor", 1.0))
        add_offset = float(_get_nc_attr(var_ord, "add_offset", 0.0))

        raw_intensity = var_ord[:].copy()
        delay_time = float(nc.variables["actual_delay_time"][:][0])
        sampling_interval = float(nc.variables["actual_sampling_interval"][:][0])
        nc.close()
    except Exception:
        import netCDF4

        nc = netCDF4.Dataset(filename, memory=file_bytes)
        # netCDF4 aplica scale_factor/add_offset sozinho por padrão. Como a escala é aplicada
        # manualmente logo abaixo, o auto-scale é desligado para não escalar duas vezes.
        nc.set_auto_maskandscale(False)
        var_ord = nc.variables["ordinate_values"]
        scale_factor = float(np.asarray(getattr(var_ord, "scale_factor", 1.0)).reshape(-1)[0])
        add_offset = float(np.asarray(getattr(var_ord, "add_offset", 0.0)).reshape(-1)[0])
        raw_intensity = np.array(var_ord[:])
        delay_time = float(np.asarray(nc.variables["actual_delay_time"][:]).reshape(-1)[0])
        sampling_interval = float(np.asarray(nc.variables["actual_sampling_interval"][:]).reshape(-1)[0])
        nc.close()

    # Aplicação do fator de escala real do detector
    intensity = (np.array(raw_intensity, dtype=np.float64) * scale_factor) + add_offset
    n_points = len(intensity)

    time_seconds = delay_time + (np.arange(n_points) * sampling_interval)
    time_minutes = time_seconds / 60.0

    return {
        "filename": filename,
        "n_points": n_points,
        "delay_time_sec": delay_time,
        "sampling_interval_sec": sampling_interval,
        "scale_factor_applied": scale_factor,
        "time": np.round(time_minutes, 4).tolist(),
        "intensity": np.round(intensity, 4).tolist(),
    }
