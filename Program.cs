using BenchmarkDotNet.Running;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BenchmarkDotNet.Configs;
using BenchmarkDotNet.Validators;
using BenchmarkDotNet.Jobs;
using BenchmarkDotNet.Environments;
using BenchmarkDotNet.Horology;
using BenchmarkDotNet.Exporters.Csv;

namespace ActivationFnBenchmarks
{
    public class Program
    {
        static void Main(string[] args)
        {
            var config = ManualConfig.Create(DefaultConfig.Instance);
            config.Add(new CsvExporter(CsvSeparator.CurrentCulture,
                new BenchmarkDotNet.Reports.SummaryStyle
                {
                    PrintUnitsInHeader = true,
                    PrintUnitsInContent = false,
                    TimeUnit = TimeUnit.Microsecond,
                    SizeUnit = BenchmarkDotNet.Columns.SizeUnit.KB
                }));

            config.Add(new Job(EnvMode.LegacyJitX64, EnvMode.Clr, RunMode.Short)
            {
                Env = { Runtime = Runtime.Clr, Platform = Platform.X64 },
                Run = { LaunchCount = 1, WarmupCount = 1, TargetCount = 1, RunStrategy = BenchmarkDotNet.Engines.RunStrategy.Throughput },
                Accuracy = { RemoveOutliers = true }
            }.WithGcAllowVeryLargeObjects(true));

            config.Add(new Job(EnvMode.RyuJitX64, EnvMode.Clr, RunMode.Short)
            {
                Env = { Runtime = Runtime.Clr, Platform = Platform.X64 },
                Run = { LaunchCount = 1, WarmupCount = 1, TargetCount = 1, RunStrategy = BenchmarkDotNet.Engines.RunStrategy.Throughput },
                Accuracy = { RemoveOutliers = true }
            }.WithGcAllowVeryLargeObjects(true));

            config.Add(BenchmarkDotNet.Loggers.ConsoleLogger.Default);
            config.Add(JitOptimizationsValidator.DontFailOnError);

            var summary = BenchmarkRunner.Run<ActivationFunctions>(config);
        }
    }
}
