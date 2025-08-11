import psutil


def get_mem(frac=0.8):
    # Determine available memory
    return int(psutil.virtual_memory().available // 1024 // 1024 * frac)
