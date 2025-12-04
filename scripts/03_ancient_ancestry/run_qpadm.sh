#!/bin/bash
#
# run_qpadm.sh - Run qpAdm analysis with multiple model configurations
#
# This script runs qpAdm with different source population models to estimate
# ancient ancestry proportions in a target sample/population.
#
# Usage:
#     bash run_qpadm.sh <qpadm_dir>
#     bash run_qpadm.sh results/qpadm/
#
# Prerequisites:
#     - ADMIXTOOLS must be installed and qpAdm must be in PATH
#     - Parameter files must be prepared using prepare_qpadm.py
#
# Models run:
#     - 2-way models (e.g., Yamnaya + Anatolia_N)
#     - 3-way models (Yamnaya + Anatolia_N + WHG)
#
# Output:
#     - {model}_output.txt - qpAdm output for each model
#     - qpadm_results_summary.txt - Summary of all results
#

set -e

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print banner
echo "========================================================"
echo "                  qpAdm Analysis Runner                  "
echo "========================================================"
echo ""

# Check arguments
if [ $# -lt 1 ]; then
    echo -e "${YELLOW}Usage: $0 <qpadm_directory>${NC}"
    echo ""
    echo "Example:"
    echo "    $0 results/qpadm/"
    echo ""
    exit 1
fi

QPADM_DIR="$1"

# Validate directory
if [ ! -d "$QPADM_DIR" ]; then
    echo -e "${RED}Error: Directory not found: $QPADM_DIR${NC}"
    exit 1
fi

# Change to qpAdm directory
cd "$QPADM_DIR"
echo -e "${BLUE}Working directory: $(pwd)${NC}"
echo ""

# Check for qpAdm executable
if ! command -v qpAdm &> /dev/null; then
    echo -e "${YELLOW}Warning: qpAdm not found in PATH${NC}"
    echo "Attempting to find qpAdm in common locations..."
    
    # Check common installation locations
    POSSIBLE_PATHS=(
        "$HOME/tools/AdmixTools/bin/qpAdm"
        "$HOME/tools/admixtools/bin/qpAdm"
        "$HOME/bin/qpAdm"
        "/usr/local/bin/qpAdm"
        "/opt/AdmixTools/bin/qpAdm"
    )
    
    QPADM_BIN=""
    for path in "${POSSIBLE_PATHS[@]}"; do
        if [ -x "$path" ]; then
            QPADM_BIN="$path"
            break
        fi
    done
    
    if [ -z "$QPADM_BIN" ]; then
        echo -e "${RED}Error: qpAdm not found. Please install ADMIXTOOLS.${NC}"
        echo ""
        echo "Installation instructions:"
        echo "  git clone https://github.com/DReichLab/AdmixTools.git"
        echo "  cd AdmixTools/src"
        echo "  make"
        echo "  # Add bin directory to PATH"
        echo ""
        exit 1
    fi
    
    echo -e "${GREEN}Found qpAdm: $QPADM_BIN${NC}"
else
    QPADM_BIN="qpAdm"
fi

# Find all parameter files
PAR_FILES=(*.par)

if [ ${#PAR_FILES[@]} -eq 0 ] || [ ! -f "${PAR_FILES[0]}" ]; then
    echo -e "${RED}Error: No parameter files (*.par) found in $QPADM_DIR${NC}"
    echo ""
    echo "Run prepare_qpadm.py first to create parameter files:"
    echo "    python scripts/03_ancient_ancestry/prepare_qpadm.py \\"
    echo "        --input data/processed/merged \\"
    echo "        --output-dir $QPADM_DIR \\"
    echo "        --sample-id SAMPLE001"
    echo ""
    exit 1
fi

echo "Found ${#PAR_FILES[@]} model configuration(s):"
for par in "${PAR_FILES[@]}"; do
    echo "  - ${par%.par}"
done
echo ""

# Initialize results summary
SUMMARY_FILE="qpadm_results_summary.txt"
echo "qpAdm Analysis Results Summary" > "$SUMMARY_FILE"
echo "==============================" >> "$SUMMARY_FILE"
echo "Date: $(date)" >> "$SUMMARY_FILE"
echo "Directory: $(pwd)" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

# Run each model
SUCCESSFUL=0
FAILED=0

for PAR_FILE in "${PAR_FILES[@]}"; do
    MODEL_NAME="${PAR_FILE%.par}"
    OUTPUT_FILE="${MODEL_NAME}_output.txt"
    
    echo -e "${BLUE}Running model: $MODEL_NAME${NC}"
    
    # Run qpAdm
    if $QPADM_BIN -p "$PAR_FILE" > "$OUTPUT_FILE" 2>&1; then
        echo -e "${GREEN}  ✓ Completed: $MODEL_NAME${NC}"
        SUCCESSFUL=$((SUCCESSFUL + 1))
        
        # Extract key results
        echo "Model: $MODEL_NAME" >> "$SUMMARY_FILE"
        echo "-" | head -c ${#MODEL_NAME} | tr '-' '-'; echo "-------" >> "$SUMMARY_FILE"
        
        # Extract coefficients and p-value
        if grep -q "best coefficients" "$OUTPUT_FILE"; then
            echo "Coefficients:" >> "$SUMMARY_FILE"
            grep -A 5 "best coefficients" "$OUTPUT_FILE" | tail -n +2 >> "$SUMMARY_FILE"
        fi
        
        if grep -q "tail prob" "$OUTPUT_FILE" || grep -q "pvalue" "$OUTPUT_FILE"; then
            echo "" >> "$SUMMARY_FILE"
            grep -E "(tail prob|pvalue|chisq|Dscore)" "$OUTPUT_FILE" >> "$SUMMARY_FILE" || true
        fi
        
        echo "" >> "$SUMMARY_FILE"
        
    else
        echo -e "${RED}  ✗ Failed: $MODEL_NAME${NC}"
        echo "    Check $OUTPUT_FILE for details"
        FAILED=$((FAILED + 1))
        
        echo "Model: $MODEL_NAME - FAILED" >> "$SUMMARY_FILE"
        echo "" >> "$SUMMARY_FILE"
    fi
    
    echo ""
done

# Print summary
echo "========================================================"
echo "                      Summary                           "
echo "========================================================"
echo -e "${GREEN}Successful: $SUCCESSFUL${NC}"
echo -e "${RED}Failed: $FAILED${NC}"
echo ""
echo -e "Results summary: ${BLUE}$QPADM_DIR/$SUMMARY_FILE${NC}"
echo ""

# Add final summary to file
echo "========================================" >> "$SUMMARY_FILE"
echo "Summary: $SUCCESSFUL successful, $FAILED failed" >> "$SUMMARY_FILE"

# Parse and display key results
echo "Key Results:"
echo "------------"

for PAR_FILE in "${PAR_FILES[@]}"; do
    MODEL_NAME="${PAR_FILE%.par}"
    OUTPUT_FILE="${MODEL_NAME}_output.txt"
    
    if [ -f "$OUTPUT_FILE" ]; then
        echo ""
        echo -e "${BLUE}$MODEL_NAME:${NC}"
        
        # Try to extract ancestry proportions
        if grep -q "best coefficients" "$OUTPUT_FILE"; then
            grep -A 10 "best coefficients" "$OUTPUT_FILE" | head -n 6
        elif grep -q "coefficients" "$OUTPUT_FILE"; then
            grep -A 5 "coefficients" "$OUTPUT_FILE" | head -n 4
        else
            echo "  (No coefficients found - check output file)"
        fi
        
        # Extract p-value if available
        if grep -q "tail prob" "$OUTPUT_FILE"; then
            grep "tail prob" "$OUTPUT_FILE" | head -n 1
        elif grep -q "pvalue" "$OUTPUT_FILE"; then
            grep "pvalue" "$OUTPUT_FILE" | head -n 1
        fi
    fi
done

echo ""
echo "========================================================"
echo "For detailed results, examine individual output files:"
for PAR_FILE in "${PAR_FILES[@]}"; do
    echo "  - ${PAR_FILE%.par}_output.txt"
done
echo "========================================================"

# Exit with appropriate code
if [ $FAILED -gt 0 ]; then
    exit 1
fi

exit 0
