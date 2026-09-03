#!/bin/bash
# Crop a PDF to its content bounding box (removes the oversized-canvas white margins) and
# re-render a matching 300-dpi PNG. Usage: crop_pdf.sh <file.pdf>
set -euo pipefail
pdf=$1; base=${pdf%.pdf}; m=6
read x0 y0 x1 y1 < <(gs -q -dBATCH -dNOPAUSE -sDEVICE=bbox "$pdf" 2>&1 \
  | awk '/^%%BoundingBox/{print $2,$3,$4,$5; exit}')
X0=$((x0 - m)); Y0=$((y0 - m)); X1=$((x1 + m)); Y1=$((y1 + m))
W=$((X1 - X0)); H=$((Y1 - Y0))
gs -q -o "${base}.crop.pdf" -sDEVICE=pdfwrite \
   -dDEVICEWIDTHPOINTS=$W -dDEVICEHEIGHTPOINTS=$H -dFIXEDMEDIA \
   -c "<</PageOffset [ $((-X0)) $((-Y0)) ]>> setpagedevice" -f "$pdf"
mv -f "${base}.crop.pdf" "$pdf"
pdftoppm -png -r 300 -singlefile "$pdf" "$base"
echo "cropped $(basename "$pdf") to ${W}x${H} pt"
