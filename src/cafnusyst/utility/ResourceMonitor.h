#pragma once

// Lightweight wall-clock/CPU-time/peak-RSS instrumentation.
//
// Usage for a one-shot step:
//   {
//     cafnusyst::ScopedResourceReport r("Configuring response_helper");
//     wu.SetResponseHelper(fclname);
//   } // prints elapsed wall/cpu time and peak RSS on scope exit
//
// Usage for a step that repeats many times (e.g. once per SRTrueInteraction),
// where per-call printing would flood stdout: accumulate and report once.
//   cafnusyst::ResourceAccumulator acc("Evaluating+saving reweights");
//   for (...) {
//     acc.Start();
//     ... do the work ...
//     acc.Stop();
//   }
//   acc.Report();

#include <chrono>
#include <cstdio>
#include <string>
#include <algorithm>
#include <sys/resource.h>

namespace cafnusyst {

struct ResourceSnapshot {
  double wall_s = 0.;    // monotonic wall-clock seconds, only meaningful as a difference
  double cpu_s = 0.;     // user+sys CPU seconds consumed by this process so far
  long maxrss_kb = 0;    // peak resident set size so far (kB on Linux); monotonically non-decreasing

  static ResourceSnapshot Now() {
    ResourceSnapshot s;
    s.wall_s = std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();

    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    s.cpu_s = (ru.ru_utime.tv_sec + ru.ru_utime.tv_usec * 1e-6) +
              (ru.ru_stime.tv_sec + ru.ru_stime.tv_usec * 1e-6);
    s.maxrss_kb = ru.ru_maxrss;
    return s;
  }
};

// RAII: reports the wall/cpu time spent and the process' peak RSS
// (as of scope exit) when it goes out of scope. Good for steps that
// run once, e.g. configuring the response_helper or building the
// global tree.
class ScopedResourceReport {
public:
  explicit ScopedResourceReport(std::string label)
      : fLabel(std::move(label)), fStart(ResourceSnapshot::Now()) {}

  ~ScopedResourceReport() { Report(); }

  void Report() {
    if (fReported) return;
    fReported = true;
    ResourceSnapshot end = ResourceSnapshot::Now();
    printf("[ResourceMonitor] %-40s wall = %9.3f s, cpu = %9.3f s, peak RSS = %9.1f MB\n",
           fLabel.c_str(), end.wall_s - fStart.wall_s, end.cpu_s - fStart.cpu_s,
           end.maxrss_kb / 1024.0);
  }

private:
  std::string fLabel;
  ResourceSnapshot fStart;
  bool fReported = false;
};

// Accumulates wall/cpu time across many short Start()/Stop() brackets
// (e.g. once per SRTrueInteraction) and reports the totals plus the
// process' peak RSS once, via Report(). Avoids flooding stdout with a
// line per event.
class ResourceAccumulator {
public:
  explicit ResourceAccumulator(std::string label = "") : fLabel(std::move(label)) {}

  void SetLabel(std::string label) { fLabel = std::move(label); }

  void Start() { fStart = ResourceSnapshot::Now(); }

  void Stop() {
    ResourceSnapshot end = ResourceSnapshot::Now();
    fWallTotal += end.wall_s - fStart.wall_s;
    fCpuTotal += end.cpu_s - fStart.cpu_s;
    fPeakRssKb = std::max(fPeakRssKb, end.maxrss_kb);
    ++fCount;
  }

  void Report() const {
    printf("[ResourceMonitor] %-40s calls = %9zu, total wall = %9.3f s (%.4f ms/call), "
           "total cpu = %9.3f s, peak RSS = %9.1f MB\n",
           fLabel.c_str(), fCount, fWallTotal, fCount ? (fWallTotal * 1000.0 / fCount) : 0.0,
           fCpuTotal, fPeakRssKb / 1024.0);
  }

  size_t Count() const { return fCount; }
  double WallTotal() const { return fWallTotal; }
  double CpuTotal() const { return fCpuTotal; }
  long PeakRssKb() const { return fPeakRssKb; }

private:
  std::string fLabel;
  ResourceSnapshot fStart;
  double fWallTotal = 0.;
  double fCpuTotal = 0.;
  long fPeakRssKb = 0;
  size_t fCount = 0;
};

} // namespace cafnusyst
