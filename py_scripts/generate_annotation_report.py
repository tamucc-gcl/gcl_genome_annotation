#!/usr/bin/env python3
"""Unified Markdown report entry and shared formatting helpers."""

def fmt(n):
    try: return f"{int(float(n)):,}"
    except (ValueError, TypeError): return str(n)

def table(header_row, data_rows):
    lines = ["| " + " | ".join(header_row) + " |",
             "|:---|" + "|".join("---:" for _ in header_row[1:]) + "|"]
    for row in data_rows:
        lines.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(lines) + "\n"

def main():
    from unified_report import main as report_main
    report_main()

if __name__ == '__main__': main()
