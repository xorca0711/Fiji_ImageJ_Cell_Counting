param(
    [string]$AuditDirectory = (Join-Path (Split-Path $PSScriptRoot -Parent) 'test_runs\current\ErrorRate_Audit_20260722')
)

$ErrorActionPreference = 'Stop'
$AuditDirectory = (Resolve-Path $AuditDirectory).Path
$workbookPath = Join-Path $AuditDirectory 'Test_Run_Error_Rate_Audit.xlsx'
$sources = [ordered]@{
    'Overall' = 'overall_error_proxies.csv'
    'Marker_Rates' = 'marker_call_rates.csv'
    'Gate_Failures' = 'morphology_gate_failure_rates.csv'
    'Call_Reasons' = 'call_reason_rates.csv'
    'Call_Status' = 'call_status_rates.csv'
    'Field_Nucleus_QC' = 'field_nucleus_qc.csv'
    'Run_Inventory' = 'run_inventory.csv'
}

function Set-Cell($cell, $value) {
    if ($null -eq $value) { $cell.Value2 = ''; return }
    $number = 0.0
    if ([double]::TryParse([string]$value, [Globalization.NumberStyles]::Float,
        [Globalization.CultureInfo]::InvariantCulture, [ref]$number)) {
        $cell.Value = $number
    } else {
        $cell.Value2 = [string]$value
    }
}

$excel = $null
$workbook = $null
try {
    $excel = New-Object -ComObject Excel.Application
    $excel.Visible = $false
    $excel.DisplayAlerts = $false
    $workbook = $excel.Workbooks.Add()
    while ($workbook.Worksheets.Count -lt $sources.Count) { $workbook.Worksheets.Add() | Out-Null }
    while ($workbook.Worksheets.Count -gt $sources.Count) { $workbook.Worksheets.Item($workbook.Worksheets.Count).Delete() }

    $sheetIndex = 1
    foreach ($entry in $sources.GetEnumerator()) {
        $sheet = $workbook.Worksheets.Item($sheetIndex++)
        $sheet.Name = $entry.Key
        $rows = @(Import-Csv -LiteralPath (Join-Path $AuditDirectory $entry.Value))
        if ($rows.Count -eq 0) { continue }
        $headers = @($rows[0].PSObject.Properties.Name)
        for ($column = 0; $column -lt $headers.Count; $column++) {
            $sheet.Cells.Item(1, $column + 1).Value2 = $headers[$column]
        }
        $rowIndex = 2
        foreach ($row in $rows) {
            for ($column = 0; $column -lt $headers.Count; $column++) {
                Set-Cell $sheet.Cells.Item($rowIndex, $column + 1) $row.($headers[$column])
            }
            $rowIndex++
        }
        $headerRange = $sheet.Range($sheet.Cells.Item(1, 1), $sheet.Cells.Item(1, $headers.Count))
        $headerRange.Font.Bold = $true
        $headerRange.Font.Color = 0xFFFFFF
        $headerRange.Interior.Color = 0x8B4513
        $headerRange.WrapText = $true
        $headerRange.AutoFilter() | Out-Null
        for ($column = 0; $column -lt $headers.Count; $column++) {
            if ($headers[$column] -match '(^percentage$|_pct|fraction)') {
                $sheet.Columns.Item($column + 1).NumberFormat = '0.00'
            }
        }
        $sheet.UsedRange.Columns.AutoFit() | Out-Null
        foreach ($column in 1..$headers.Count) {
            if ($sheet.Columns.Item($column).ColumnWidth -gt 55) { $sheet.Columns.Item($column).ColumnWidth = 55 }
        }
        $sheet.UsedRange.Rows.AutoFit() | Out-Null
    }

    $workbook.SaveAs($workbookPath, 51)
    $workbook.Close($true)
    $excel.Quit()
} finally {
    if ($workbook -ne $null) { [void][Runtime.InteropServices.Marshal]::FinalReleaseComObject($workbook) }
    if ($excel -ne $null) { [void][Runtime.InteropServices.Marshal]::FinalReleaseComObject($excel) }
    [GC]::Collect(); [GC]::WaitForPendingFinalizers()
}

Write-Output $workbookPath
