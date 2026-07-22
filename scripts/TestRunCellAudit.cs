using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

public static class TestRunCellAudit
{
    private sealed class GateCount { public long Pass, Fail, Missing; }
    private sealed class Metric
    {
        public long Opportunities, Evaluable, Positive, Negative, Indeterminate;
        public long RawPositive, RawNegative, ConcordantPositive, ConcordantNegative;
        public long RawPositiveFinalNegative, RawNegativeFinalPositive;
        public long RawPositiveIndeterminate, RawNegativeIndeterminate;
        public readonly Dictionary<string, long> Statuses = new Dictionary<string, long>();
        public readonly Dictionary<string, long> Reasons = new Dictionary<string, long>();
        public readonly Dictionary<string, GateCount> Gates = new Dictionary<string, GateCount>();
    }

    private static readonly string[] GateNames = {
        "fraction_pass", "connected_pattern_pass", "ownership_clear",
        "projection_valid", "compartment_pass", "enrichment_pass"
    };

    private static void Increment(Dictionary<string, long> table, string key, long amount = 1)
    {
        long old; table.TryGetValue(key, out old); table[key] = old + amount;
    }

    private static Metric GetMetric(Dictionary<string, Metric> table, string key)
    {
        Metric value;
        if (!table.TryGetValue(key, out value)) { value = new Metric(); table[key] = value; }
        return value;
    }

    private static void Add(Metric target, Metric source)
    {
        target.Opportunities += source.Opportunities; target.Evaluable += source.Evaluable;
        target.Positive += source.Positive; target.Negative += source.Negative;
        target.Indeterminate += source.Indeterminate; target.RawPositive += source.RawPositive;
        target.RawNegative += source.RawNegative; target.ConcordantPositive += source.ConcordantPositive;
        target.ConcordantNegative += source.ConcordantNegative;
        target.RawPositiveFinalNegative += source.RawPositiveFinalNegative;
        target.RawNegativeFinalPositive += source.RawNegativeFinalPositive;
        target.RawPositiveIndeterminate += source.RawPositiveIndeterminate;
        target.RawNegativeIndeterminate += source.RawNegativeIndeterminate;
        foreach (var pair in source.Statuses) Increment(target.Statuses, pair.Key, pair.Value);
        foreach (var pair in source.Reasons) Increment(target.Reasons, pair.Key, pair.Value);
        foreach (var pair in source.Gates)
        {
            GateCount gate;
            if (!target.Gates.TryGetValue(pair.Key, out gate)) { gate = new GateCount(); target.Gates[pair.Key] = gate; }
            gate.Pass += pair.Value.Pass; gate.Fail += pair.Value.Fail; gate.Missing += pair.Value.Missing;
        }
    }

    private static string Csv(object value)
    {
        if (value == null) return "";
        string text = Convert.ToString(value, CultureInfo.InvariantCulture) ?? "";
        if (text.IndexOfAny(new[] { ',', '"', '\r', '\n' }) >= 0) return "\"" + text.Replace("\"", "\"\"") + "\"";
        return text;
    }

    private static string Pct(long numerator, long denominator)
    {
        return denominator == 0 ? "" : (100.0 * numerator / denominator).ToString("R", CultureInfo.InvariantCulture);
    }

    private static void Update(Metric m, string raw, string finalCall, string status, string reason,
                               string[] values, Dictionary<string, int> gateIndexes)
    {
        m.Opportunities++;
        if (raw == "1") m.RawPositive++; else m.RawNegative++;
        Increment(m.Statuses, String.IsNullOrWhiteSpace(status) ? "missing" : status);
        if (finalCall == "1")
        {
            m.Evaluable++; m.Positive++;
            if (raw == "1") m.ConcordantPositive++; else m.RawNegativeFinalPositive++;
        }
        else if (finalCall == "0")
        {
            m.Evaluable++; m.Negative++;
            if (raw == "0") m.ConcordantNegative++; else m.RawPositiveFinalNegative++;
        }
        else
        {
            m.Indeterminate++;
            if (raw == "1") m.RawPositiveIndeterminate++; else m.RawNegativeIndeterminate++;
        }
        if (!String.IsNullOrWhiteSpace(reason))
            foreach (string token in reason.Split(';')) if (!String.IsNullOrWhiteSpace(token)) Increment(m.Reasons, token.Trim());
        foreach (var pair in gateIndexes)
        {
            GateCount gate;
            if (!m.Gates.TryGetValue(pair.Key, out gate)) { gate = new GateCount(); m.Gates[pair.Key] = gate; }
            string value = pair.Value < values.Length ? values[pair.Value] : "";
            if (value == "1") gate.Pass++; else if (value == "0") gate.Fail++; else gate.Missing++;
        }
    }

    public static void Run(string root, string output)
    {
        Directory.CreateDirectory(output);
        var sets = new Dictionary<string, Dictionary<string, Metric>> {
            { "Canonical45", new Dictionary<string, Metric>() },
            { "AllAttempts47", new Dictionary<string, Metric>() }
        };

        foreach (string file in Directory.EnumerateFiles(root, "*_cells.csv", SearchOption.AllDirectories))
        {
            bool canonical = file.IndexOf("FinalPilot_", StringComparison.OrdinalIgnoreCase) < 0;
            using (var reader = new StreamReader(file, Encoding.UTF8, true, 1 << 16))
            {
                string headerLine = reader.ReadLine();
                if (headerLine == null) continue;
                string[] headers = headerLine.Split(new[] { ',' }, StringSplitOptions.None);
                var index = headers.Select((name, i) => new { name, i }).ToDictionary(x => x.name, x => x.i);
                int panelIndex = index["panel"];
                var markerInfo = new Dictionary<string, Tuple<int, int, int, int, Dictionary<string, int>>>();
                foreach (string marker in new[] { "CC10", "T1A", "tdTOM", "mRAGE", "AcTub" })
                {
                    int statusIndex;
                    if (!index.TryGetValue(marker + "_call_status", out statusIndex)) continue;
                    var gateIndexes = new Dictionary<string, int>();
                    foreach (string gate in GateNames)
                    {
                        int gateIndex;
                        if (index.TryGetValue(marker + "_" + gate, out gateIndex)) gateIndexes[gate] = gateIndex;
                    }
                    markerInfo[marker] = Tuple.Create(index[marker + "_pos"], index[marker + "_final_call"],
                                                      statusIndex, index[marker + "_call_reason"], gateIndexes);
                }

                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    if (line.Length == 0) continue;
                    string[] values = line.Split(new[] { ',' }, StringSplitOptions.None);
                    string panel = values[panelIndex];
                    foreach (var markerPair in markerInfo)
                    {
                        var mi = markerPair.Value;
                        foreach (string setName in canonical ? new[] { "Canonical45", "AllAttempts47" } : new[] { "AllAttempts47" })
                        {
                            Metric m = GetMetric(sets[setName], panel + "|" + markerPair.Key);
                            Update(m, values[mi.Item1], values[mi.Item2], values[mi.Item3], values[mi.Item4], values, mi.Item5);
                        }
                    }
                }
            }
        }

        foreach (var setPair in sets)
        {
            var aggregateSpecs = new Dictionary<string, string[]> {
                { "E|ALL_ASSIGNABLE", new[] { "E|CC10", "E|tdTOM" } },
                { "R|ALL_ASSIGNABLE", new[] { "R|T1A", "R|tdTOM", "R|mRAGE" } },
                { "ALL|tdTOM", new[] { "E|tdTOM", "R|tdTOM" } },
                { "ALL|ALL_ASSIGNABLE", new[] { "E|CC10", "E|tdTOM", "R|T1A", "R|tdTOM", "R|mRAGE" } }
            };
            foreach (var spec in aggregateSpecs)
            {
                var aggregate = new Metric();
                foreach (string sourceKey in spec.Value) { Metric source; if (setPair.Value.TryGetValue(sourceKey, out source)) Add(aggregate, source); }
                setPair.Value[spec.Key] = aggregate;
            }
        }

        using (var marker = new StreamWriter(Path.Combine(output, "marker_call_rates.csv"), false, new UTF8Encoding(true)))
        using (var status = new StreamWriter(Path.Combine(output, "call_status_rates.csv"), false, new UTF8Encoding(true)))
        using (var reason = new StreamWriter(Path.Combine(output, "call_reason_rates.csv"), false, new UTF8Encoding(true)))
        using (var gate = new StreamWriter(Path.Combine(output, "morphology_gate_failure_rates.csv"), false, new UTF8Encoding(true)))
        {
            marker.WriteLine("dataset,panel,marker,opportunities,evaluable,evaluable_pct,indeterminate,indeterminate_pct,true_positive,true_positive_pct_of_evaluable,true_negative,true_negative_pct_of_evaluable,raw_intensity_positive,raw_intensity_positive_pct,concordant,concordance_pct_of_evaluable,raw_positive_final_negative,raw_negative_final_positive,intensity_morphology_discordant,discordance_pct_of_evaluable,raw_positive_indeterminate,raw_negative_indeterminate,review_burden_proxy,review_burden_pct_of_opportunities,confirmatory_calls,confirmatory_pct_of_evaluable,exploratory_calls,exploratory_pct_of_evaluable");
            status.WriteLine("dataset,panel,marker,status,count,pct_of_opportunities");
            reason.WriteLine("dataset,panel,marker,reason,occurrences,pct_of_opportunities");
            gate.WriteLine("dataset,panel,marker,gate,pass,fail,missing,fail_pct_of_known");
            foreach (var setPair in sets.OrderBy(x => x.Key))
            foreach (var metricPair in setPair.Value.OrderBy(x => x.Key))
            {
                string[] parts = metricPair.Key.Split('|'); Metric m = metricPair.Value;
                long concordant = m.ConcordantPositive + m.ConcordantNegative;
                long discordant = m.RawPositiveFinalNegative + m.RawNegativeFinalPositive;
                long confirmatory = Get(m.Statuses, "positive") + Get(m.Statuses, "negative");
                long exploratory = Get(m.Statuses, "exploratory_positive") + Get(m.Statuses, "exploratory_negative");
                object[] values = { setPair.Key, parts[0], parts[1], m.Opportunities, m.Evaluable, Pct(m.Evaluable,m.Opportunities),
                    m.Indeterminate, Pct(m.Indeterminate,m.Opportunities), m.Positive, Pct(m.Positive,m.Evaluable),
                    m.Negative, Pct(m.Negative,m.Evaluable), m.RawPositive, Pct(m.RawPositive,m.Opportunities),
                    concordant, Pct(concordant,m.Evaluable), m.RawPositiveFinalNegative, m.RawNegativeFinalPositive,
                    discordant, Pct(discordant,m.Evaluable), m.RawPositiveIndeterminate, m.RawNegativeIndeterminate,
                    m.Indeterminate + discordant, Pct(m.Indeterminate+discordant,m.Opportunities),
                    confirmatory, Pct(confirmatory,m.Evaluable), exploratory, Pct(exploratory,m.Evaluable) };
                marker.WriteLine(String.Join(",", values.Select(Csv)));
                foreach (var pair in m.Statuses.OrderBy(x => x.Key))
                    status.WriteLine(String.Join(",", new object[] { setPair.Key,parts[0],parts[1],pair.Key,pair.Value,Pct(pair.Value,m.Opportunities) }.Select(Csv)));
                foreach (var pair in m.Reasons.OrderBy(x => x.Key))
                    reason.WriteLine(String.Join(",", new object[] { setPair.Key,parts[0],parts[1],pair.Key,pair.Value,Pct(pair.Value,m.Opportunities) }.Select(Csv)));
                foreach (var pair in m.Gates.OrderBy(x => x.Key))
                    gate.WriteLine(String.Join(",", new object[] { setPair.Key,parts[0],parts[1],pair.Key,pair.Value.Pass,pair.Value.Fail,pair.Value.Missing,Pct(pair.Value.Fail,pair.Value.Pass+pair.Value.Fail) }.Select(Csv)));
            }
        }
    }

    private static long Get(Dictionary<string, long> table, string key) { long value; return table.TryGetValue(key, out value) ? value : 0; }
}
