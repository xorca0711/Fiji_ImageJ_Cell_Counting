param(
    [string]$TestRunsRoot = (Join-Path (Split-Path $PSScriptRoot -Parent) 'test_runs'),
    [string]$OutputDirectory = (Join-Path (Split-Path $PSScriptRoot -Parent) 'test_runs\current\ErrorRate_Audit_20260722')
)

$ErrorActionPreference = 'Stop'
$TestRunsRoot = (Resolve-Path $TestRunsRoot).Path
New-Item -ItemType Directory -Path $OutputDirectory -Force | Out-Null
$OutputDirectory = (Resolve-Path $OutputDirectory).Path

if (-not ('TestRunCellAudit' -as [type])) {
    Add-Type -Path (Join-Path $PSScriptRoot 'TestRunCellAudit.cs')
}
[TestRunCellAudit]::Run($TestRunsRoot, $OutputDirectory)

function Number([object]$value) {
    if ($null -eq $value -or [string]::IsNullOrWhiteSpace([string]$value)) { return 0.0 }
    return [double]::Parse([string]$value, [Globalization.CultureInfo]::InvariantCulture)
}
function Percent([double]$numerator, [double]$denominator) {
    if ($denominator -le 0) { return $null }
    return 100.0 * $numerator / $denominator
}

$summary = @()
$inventory = @()
foreach ($file in Get-ChildItem -LiteralPath $TestRunsRoot -Recurse -Filter 'run_summary.csv') {
    $relative = $file.FullName.Substring($TestRunsRoot.Length + 1)
    $canonical = $relative -notlike '*FinalPilot_*'
    $rows = @(Import-Csv -LiteralPath $file.FullName)
    $inventory += [pscustomobject]@{
        run_summary = $relative; rows = $rows.Count
        dataset = if ($canonical) { 'Canonical45' } else { 'PilotDuplicate' }
    }
    foreach ($row in $rows) {
        $summary += [pscustomobject]@{
            dataset = if ($canonical) { 'Canonical45' } else { 'PilotDuplicate' }
            canonical = $canonical; run_summary = $relative; image = $row.image
            panel = $row.panel; condition = $row.condition; section_id = $row.section_id
            accepted_nuclei = [long](Number $row.n_nuclei)
            rejected_candidates = [long](Number $row.n_rejected_nucleus_candidates)
            rejected_below_min_area = [long](Number $row.n_rejected_below_min_area)
            rejected_at_image_edge = [long](Number $row.n_rejected_at_image_edge)
            rejected_by_particle_filter = [long](Number $row.n_rejected_by_particle_filter)
        }
    }
}

$fieldRows = foreach ($row in $summary) {
    $candidates = $row.accepted_nuclei + $row.rejected_candidates
    $flags = @()
    if ($row.accepted_nuclei -lt 1500) { $flags += 'low_accepted_nuclei' }
    if ($row.rejected_candidates / [Math]::Max(1.0, $row.accepted_nuclei) -ge 2.0) { $flags += 'extreme_rejected_to_accepted_ratio' }
    [pscustomobject]@{
        dataset = $row.dataset; canonical = $row.canonical; run_summary = $row.run_summary
        image = $row.image; panel = $row.panel; condition = $row.condition; section_id = $row.section_id
        accepted_nuclei = $row.accepted_nuclei; rejected_candidates = $row.rejected_candidates
        total_candidates = $candidates; acceptance_pct = Percent $row.accepted_nuclei $candidates
        rejection_pct = Percent $row.rejected_candidates $candidates
        rejected_below_min_area = $row.rejected_below_min_area
        below_min_pct_of_rejected = Percent $row.rejected_below_min_area $row.rejected_candidates
        rejected_at_image_edge = $row.rejected_at_image_edge
        edge_pct_of_rejected = Percent $row.rejected_at_image_edge $row.rejected_candidates
        rejected_by_particle_filter = $row.rejected_by_particle_filter
        qc_flag = ($flags -join ';')
    }
}
$fieldRows | Export-Csv -LiteralPath (Join-Path $OutputDirectory 'field_nucleus_qc.csv') -NoTypeInformation -Encoding UTF8
$inventory | Export-Csv -LiteralPath (Join-Path $OutputDirectory 'run_inventory.csv') -NoTypeInformation -Encoding UTF8

$markerRows = @(Import-Csv -LiteralPath (Join-Path $OutputDirectory 'marker_call_rates.csv'))
$overview = @()
foreach ($setName in @('Canonical45', 'AllAttempts47')) {
    $fieldSet = if ($setName -eq 'Canonical45') { @($summary | Where-Object canonical) } else { @($summary) }
    $accepted = [long](($fieldSet | Measure-Object accepted_nuclei -Sum).Sum)
    $rejected = [long](($fieldSet | Measure-Object rejected_candidates -Sum).Sum)
    $below = [long](($fieldSet | Measure-Object rejected_below_min_area -Sum).Sum)
    $edge = [long](($fieldSet | Measure-Object rejected_at_image_edge -Sum).Sum)
    $particle = [long](($fieldSet | Measure-Object rejected_by_particle_filter -Sum).Sum)
    $candidateTotal = $accepted + $rejected
    $m = $markerRows | Where-Object { $_.dataset -eq $setName -and $_.panel -eq 'ALL' -and $_.marker -eq 'ALL_ASSIGNABLE' }
    $opportunities = [long](Number $m.opportunities)
    $evaluable = [long](Number $m.evaluable)
    $indeterminate = [long](Number $m.indeterminate)
    $positive = [long](Number $m.true_positive)
    $negative = [long](Number $m.true_negative)
    $concordant = [long](Number $m.concordant)
    $discordant = [long](Number $m.intensity_morphology_discordant)
    $rawPosFinalNeg = [long](Number $m.raw_positive_final_negative)
    $rawNegFinalPos = [long](Number $m.raw_negative_final_positive)
    $confirmatory = [long](Number $m.confirmatory_calls)
    $exploratory = [long](Number $m.exploratory_calls)
    $specs = @(
        @('Inventory', 'Fields analyzed', $fieldSet.Count, $fieldSet.Count, 'Field count'),
        @('Nucleus segmentation', 'Accepted nucleus candidates', $accepted, $candidateTotal, 'Accepted / all DAPI-derived candidates'),
        @('Nucleus segmentation', 'Rejected nucleus candidates', $rejected, $candidateTotal, 'QC loss proxy; not a validated biological error rate'),
        @('Nucleus rejection reason', 'Below minimum area', $below, $rejected, 'Fraction of rejected candidates'),
        @('Nucleus rejection reason', 'At image edge', $edge, $rejected, 'Fraction of rejected candidates'),
        @('Nucleus rejection reason', 'Particle filter', $particle, $rejected, 'Fraction of rejected candidates'),
        @('Assignable marker calls', 'Evaluable calls', $evaluable, $opportunities, 'Morphology decision was possible'),
        @('Assignable marker calls', 'Indeterminate calls', $indeterminate, $opportunities, 'Non-evaluable; main marker-call QC loss proxy'),
        @('Assignable marker calls', 'True positive calls', $positive, $evaluable, 'Morphology-positive / evaluable'),
        @('Assignable marker calls', 'True negative calls', $negative, $evaluable, 'Morphology-negative / evaluable'),
        @('Intensity vs morphology', 'Concordant calls', $concordant, $evaluable, 'Raw mean intensity and morphology final call agree'),
        @('Intensity vs morphology', 'Discordant calls', $discordant, $evaluable, 'Sensitivity-analysis proxy; morphology remains authoritative'),
        @('Intensity vs morphology', 'Raw positive changed to final negative', $rawPosFinalNeg, $evaluable, 'Morphology rejected an intensity-positive call'),
        @('Intensity vs morphology', 'Raw negative changed to final positive', $rawNegFinalPos, $evaluable, 'Morphology rescued localized signal missed by mean intensity'),
        @('Review burden proxy', 'Indeterminate or intensity-morphology discordant', ($indeterminate + $discordant), $opportunities, 'Non-overlapping calls requiring exclusion or method review'),
        @('Threshold status', 'Confirmatory fixed-threshold calls', $confirmatory, $evaluable, 'Final-study approval requires control-derived fixed thresholds'),
        @('Threshold status', 'Exploratory adaptive-threshold calls', $exploratory, $evaluable, 'Adaptive Otsu test-run calls')
    )
    foreach ($spec in $specs) {
        $overview += [pscustomobject]@{
            dataset = $setName; domain = $spec[0]; metric = $spec[1]
            numerator = $spec[2]; denominator = $spec[3]; percentage = Percent $spec[2] $spec[3]
            interpretation = $spec[4]
        }
    }
}
$overview | Export-Csv -LiteralPath (Join-Path $OutputDirectory 'overall_error_proxies.csv') -NoTypeInformation -Encoding UTF8

Write-Output ('Analyzed cell tables: ' + @(Get-ChildItem -LiteralPath $TestRunsRoot -Recurse -Filter '*_cells.csv').Count)
Write-Output ('Canonical fields: ' + @($summary | Where-Object canonical).Count)
Write-Output ('All attempts: ' + $summary.Count)
Write-Output $OutputDirectory
