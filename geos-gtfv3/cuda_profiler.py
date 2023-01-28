import time

# Conditional cupy import for non-GPU machines
try:
    import cupy as cp
except ModuleNotFoundError:
    cp = None

# Run a deviceSynchronize() to check
# that the GPU is present and ready to run
if cp is not None:
    try:
        cp.cuda.runtime.deviceSynchronize()
        GPU_AVAILABLE = True
    except cp.cuda.runtime.CUDARuntimeError:
        GPU_AVAILABLE = False
else:
    GPU_AVAILABLE = False


class CUDAProfiler:
    """Leverages NVTX & NSYS to profile CUDA kernels."""

    def __init__(self, label: str) -> None:
        self.label = label
        self._start_time = 0

    def __enter__(self):
        if cp is not None:
            cp.cuda.runtime.deviceSynchronize()
            cp.cuda.nvtx.RangePush(self.label)
        self._start_time = time.perf_counter()

    def __exit__(self, _type, _val, _traceback):
        if cp is not None:
            cp.cuda.runtime.deviceSynchronize()
            cp.cuda.nvtx.RangePop()
        t = time.perf_counter() - self._start_time
        print(f"{self.label} CPU time: {t}s")

    @classmethod
    def sync_device():
        if GPU_AVAILABLE:
            cp.cuda.runtime.deviceSynchronize()

    @classmethod
    def start_cuda_profiler(cls):
        if GPU_AVAILABLE:
            cp.cuda.profiler.start()

    @classmethod
    def stop_cuda_profiler(cls):
        if GPU_AVAILABLE:
            cp.cuda.profiler.stop()

    @classmethod
    def mark_cuda_profiler(cls, message: str):
        if GPU_AVAILABLE:
            cp.cuda.nvtx.Mark(message)
