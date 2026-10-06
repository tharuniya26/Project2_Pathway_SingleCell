import zipfile
from pathlib import Path

zip_path = "data/raw/Rmarkdown.zip"
output_path = Path("documents/pseudobulk_conversion_function.txt")

with zipfile.ZipFile(zip_path) as z:
    files = [
        name for name in z.namelist()
        if name.lower().endswith((".r", ".rmd"))
    ]

    for name in files:
        lines = z.read(name).decode(
            "utf-8",
            errors="ignore"
        ).splitlines()

        for i, line in enumerate(lines):
            if (
                "pseudobulk_conversion" in line.lower()
                and "function" in line.lower()
            ):
                start = max(0, i - 5)
                end = min(len(lines), i + 100)

                excerpt = [
                    f"SOURCE FILE: {name}",
                    f"FUNCTION START LINE: {i + 1}",
                    "=" * 70,
                ]

                for j in range(start, end):
                    excerpt.append(f"{j + 1}: {lines[j]}")

                output_path.write_text(
                    "\n".join(excerpt),
                    encoding="utf-8"
                )

                print("Saved function to:")
                print(output_path)

                raise SystemExit