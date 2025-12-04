#!/bin/bash
#
# install_tools.sh - Install bioinformatics tools for ancestry analysis
#
# This script installs the following tools:
#   - PLINK 1.9 (download pre-compiled binary)
#   - ADMIXTURE (download from website)
#   - EIGENSOFT (compile from source or download)
#   - ADMIXTOOLS (compile from source)
#   - Optional: RFMix, TreeMix
#
# Usage:
#     bash install_tools.sh [--install-dir ~/tools]
#
# Requirements:
#     - wget or curl for downloading
#     - gcc, make, gsl, openblas for compilation
#     - Ubuntu/Debian: apt-get install build-essential libgsl-dev libopenblas-dev
#     - macOS: brew install gsl openblas
#

set -e

# Default installation directory
INSTALL_DIR="${HOME}/tools"

# Tool versions
PLINK_VERSION="1.9"
ADMIXTURE_VERSION="1.3.0"
EIGENSOFT_VERSION="8.0.0"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print banner
print_banner() {
    echo "========================================================"
    echo "      Ancestry Analysis Pipeline - Tool Installer       "
    echo "========================================================"
    echo ""
}

# Print usage
print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --install-dir DIR    Installation directory (default: ~/tools)"
    echo "  --plink              Install PLINK only"
    echo "  --admixture          Install ADMIXTURE only"
    echo "  --eigensoft          Install EIGENSOFT only"
    echo "  --admixtools         Install ADMIXTOOLS only"
    echo "  --rfmix              Install RFMix"
    echo "  --treemix            Install TreeMix"
    echo "  --all                Install all tools (default)"
    echo "  --help               Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 --install-dir ~/tools --all"
    echo ""
}

# Parse command line arguments
INSTALL_PLINK=false
INSTALL_ADMIXTURE=false
INSTALL_EIGENSOFT=false
INSTALL_ADMIXTOOLS=false
INSTALL_RFMIX=false
INSTALL_TREEMIX=false
INSTALL_ALL=true

while [[ $# -gt 0 ]]; do
    case $1 in
        --install-dir)
            INSTALL_DIR="$2"
            shift 2
            ;;
        --plink)
            INSTALL_PLINK=true
            INSTALL_ALL=false
            shift
            ;;
        --admixture)
            INSTALL_ADMIXTURE=true
            INSTALL_ALL=false
            shift
            ;;
        --eigensoft)
            INSTALL_EIGENSOFT=true
            INSTALL_ALL=false
            shift
            ;;
        --admixtools)
            INSTALL_ADMIXTOOLS=true
            INSTALL_ALL=false
            shift
            ;;
        --rfmix)
            INSTALL_RFMIX=true
            INSTALL_ALL=false
            shift
            ;;
        --treemix)
            INSTALL_TREEMIX=true
            INSTALL_ALL=false
            shift
            ;;
        --all)
            INSTALL_ALL=true
            shift
            ;;
        --help)
            print_usage
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            print_usage
            exit 1
            ;;
    esac
done

# If --all, set all flags
if $INSTALL_ALL; then
    INSTALL_PLINK=true
    INSTALL_ADMIXTURE=true
    INSTALL_EIGENSOFT=true
    INSTALL_ADMIXTOOLS=true
fi

# Detect operating system
detect_os() {
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        OS="linux"
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        OS="macos"
    else
        echo -e "${RED}Unsupported operating system: $OSTYPE${NC}"
        exit 1
    fi
    echo -e "${BLUE}Detected OS: $OS${NC}"
}

# Check for required tools
check_prerequisites() {
    echo -e "${BLUE}Checking prerequisites...${NC}"
    
    local missing=()
    
    # Check for download tools
    if ! command -v wget &> /dev/null && ! command -v curl &> /dev/null; then
        missing+=("wget or curl")
    fi
    
    # Check for unzip
    if ! command -v unzip &> /dev/null; then
        missing+=("unzip")
    fi
    
    # Check for compilation tools (if needed)
    if $INSTALL_EIGENSOFT || $INSTALL_ADMIXTOOLS; then
        if ! command -v gcc &> /dev/null; then
            missing+=("gcc")
        fi
        if ! command -v make &> /dev/null; then
            missing+=("make")
        fi
    fi
    
    if [ ${#missing[@]} -gt 0 ]; then
        echo -e "${YELLOW}Missing prerequisites: ${missing[*]}${NC}"
        echo ""
        echo "Install them with:"
        if [[ "$OS" == "linux" ]]; then
            echo "  sudo apt-get update"
            echo "  sudo apt-get install -y build-essential wget unzip libgsl-dev libopenblas-dev"
        elif [[ "$OS" == "macos" ]]; then
            echo "  brew install wget gsl openblas"
        fi
        echo ""
        return 1
    fi
    
    echo -e "${GREEN}All prerequisites satisfied${NC}"
    return 0
}

# Download file
download_file() {
    local url="$1"
    local output="$2"
    
    if command -v wget &> /dev/null; then
        wget -q -O "$output" "$url"
    elif command -v curl &> /dev/null; then
        curl -sL -o "$output" "$url"
    else
        echo -e "${RED}Neither wget nor curl available${NC}"
        return 1
    fi
}

# Install PLINK 1.9
install_plink() {
    echo ""
    echo -e "${BLUE}Installing PLINK 1.9...${NC}"
    
    local plink_dir="${INSTALL_DIR}/plink"
    mkdir -p "$plink_dir"
    
    local download_url
    if [[ "$OS" == "linux" ]]; then
        download_url="https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip"
    elif [[ "$OS" == "macos" ]]; then
        download_url="https://s3.amazonaws.com/plink1-assets/plink_mac_20231211.zip"
    fi
    
    local temp_zip="${plink_dir}/plink.zip"
    
    echo "  Downloading PLINK..."
    if download_file "$download_url" "$temp_zip"; then
        echo "  Extracting..."
        unzip -q -o "$temp_zip" -d "$plink_dir"
        rm "$temp_zip"
        chmod +x "${plink_dir}/plink"
        
        echo -e "${GREEN}  ✓ PLINK installed to: ${plink_dir}/plink${NC}"
    else
        echo -e "${RED}  ✗ Failed to download PLINK${NC}"
        return 1
    fi
}

# Install ADMIXTURE
install_admixture() {
    echo ""
    echo -e "${BLUE}Installing ADMIXTURE ${ADMIXTURE_VERSION}...${NC}"
    
    local admixture_dir="${INSTALL_DIR}/admixture"
    mkdir -p "$admixture_dir"
    
    local download_url
    if [[ "$OS" == "linux" ]]; then
        download_url="https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz"
    elif [[ "$OS" == "macos" ]]; then
        download_url="https://dalexander.github.io/admixture/binaries/admixture_macosx-1.3.0.tar.gz"
    fi
    
    local temp_tar="${admixture_dir}/admixture.tar.gz"
    
    echo "  Downloading ADMIXTURE..."
    if download_file "$download_url" "$temp_tar"; then
        echo "  Extracting..."
        tar -xzf "$temp_tar" -C "$admixture_dir" --strip-components=1
        rm "$temp_tar"
        chmod +x "${admixture_dir}/admixture"
        
        echo -e "${GREEN}  ✓ ADMIXTURE installed to: ${admixture_dir}/admixture${NC}"
    else
        echo -e "${RED}  ✗ Failed to download ADMIXTURE${NC}"
        return 1
    fi
}

# Install EIGENSOFT
install_eigensoft() {
    echo ""
    echo -e "${BLUE}Installing EIGENSOFT...${NC}"
    
    local eigensoft_dir="${INSTALL_DIR}/eigensoft"
    mkdir -p "$eigensoft_dir"
    
    # Clone from GitHub
    echo "  Cloning EIGENSOFT repository..."
    if command -v git &> /dev/null; then
        local temp_dir="${eigensoft_dir}/EIG_source"
        
        if [ -d "$temp_dir" ]; then
            rm -rf "$temp_dir"
        fi
        
        git clone --quiet https://github.com/DReichLab/EIG.git "$temp_dir"
        
        echo "  Compiling EIGENSOFT..."
        cd "$temp_dir/src"
        
        # Modify CFLAGS if needed for macOS
        if [[ "$OS" == "macos" ]]; then
            # Add OpenBLAS path if using brew
            if command -v brew &> /dev/null; then
                export LDFLAGS="-L$(brew --prefix openblas)/lib"
                export CFLAGS="-I$(brew --prefix openblas)/include"
            fi
        fi
        
        if make -j2 2>/dev/null; then
            # Copy binaries
            cp smartpca convertf "$eigensoft_dir/"
            chmod +x "${eigensoft_dir}/smartpca" "${eigensoft_dir}/convertf"
            
            # Clean up
            cd - > /dev/null
            rm -rf "$temp_dir"
            
            echo -e "${GREEN}  ✓ EIGENSOFT installed to: ${eigensoft_dir}${NC}"
        else
            echo -e "${YELLOW}  Compilation failed. Trying pre-built binaries...${NC}"
            cd - > /dev/null
            rm -rf "$temp_dir"
            
            # Fall back to alternative method
            echo -e "${YELLOW}  Note: EIGENSOFT may need to be installed manually${NC}"
            echo "  See: https://github.com/DReichLab/EIG"
        fi
    else
        echo -e "${YELLOW}  git not available. Please install EIGENSOFT manually.${NC}"
        echo "  See: https://github.com/DReichLab/EIG"
    fi
}

# Install ADMIXTOOLS
install_admixtools() {
    echo ""
    echo -e "${BLUE}Installing ADMIXTOOLS...${NC}"
    
    local admixtools_dir="${INSTALL_DIR}/admixtools"
    mkdir -p "$admixtools_dir"
    
    # Clone from GitHub
    echo "  Cloning ADMIXTOOLS repository..."
    if command -v git &> /dev/null; then
        local temp_dir="${admixtools_dir}/AdmixTools_source"
        
        if [ -d "$temp_dir" ]; then
            rm -rf "$temp_dir"
        fi
        
        git clone --quiet https://github.com/DReichLab/AdmixTools.git "$temp_dir"
        
        echo "  Compiling ADMIXTOOLS..."
        cd "$temp_dir/src"
        
        # Set up environment for compilation
        if [[ "$OS" == "macos" ]]; then
            if command -v brew &> /dev/null; then
                export LDFLAGS="-L$(brew --prefix openblas)/lib -L$(brew --prefix gsl)/lib"
                export CFLAGS="-I$(brew --prefix openblas)/include -I$(brew --prefix gsl)/include"
            fi
        fi
        
        if make -j2 2>/dev/null; then
            # Copy binaries
            mkdir -p "${admixtools_dir}/bin"
            cp qpAdm qpDstat qpGraph qpF4ratio rolloff "${admixtools_dir}/bin/" 2>/dev/null || true
            chmod +x "${admixtools_dir}/bin/"* 2>/dev/null || true
            
            # Clean up
            cd - > /dev/null
            rm -rf "$temp_dir"
            
            echo -e "${GREEN}  ✓ ADMIXTOOLS installed to: ${admixtools_dir}/bin${NC}"
        else
            echo -e "${YELLOW}  Compilation failed.${NC}"
            cd - > /dev/null
            rm -rf "$temp_dir"
            
            echo -e "${YELLOW}  Note: ADMIXTOOLS may need to be installed manually${NC}"
            echo "  See: https://github.com/DReichLab/AdmixTools"
        fi
    else
        echo -e "${YELLOW}  git not available. Please install ADMIXTOOLS manually.${NC}"
    fi
}

# Install RFMix
install_rfmix() {
    echo ""
    echo -e "${BLUE}Installing RFMix...${NC}"
    
    local rfmix_dir="${INSTALL_DIR}/rfmix"
    mkdir -p "$rfmix_dir"
    
    if command -v git &> /dev/null; then
        echo "  Cloning RFMix repository..."
        
        if [ -d "${rfmix_dir}/RFMix_source" ]; then
            rm -rf "${rfmix_dir}/RFMix_source"
        fi
        
        git clone --quiet https://github.com/slowkoni/rfmix.git "${rfmix_dir}/RFMix_source"
        
        cd "${rfmix_dir}/RFMix_source"
        
        # Build
        if autoreconf --force --install 2>/dev/null && ./configure && make -j2; then
            cp rfmix "$rfmix_dir/"
            chmod +x "${rfmix_dir}/rfmix"
            cd - > /dev/null
            rm -rf "${rfmix_dir}/RFMix_source"
            
            echo -e "${GREEN}  ✓ RFMix installed to: ${rfmix_dir}/rfmix${NC}"
        else
            echo -e "${YELLOW}  RFMix compilation failed${NC}"
            cd - > /dev/null
        fi
    else
        echo -e "${YELLOW}  git not available. Please install RFMix manually.${NC}"
    fi
}

# Install TreeMix
install_treemix() {
    echo ""
    echo -e "${BLUE}Installing TreeMix...${NC}"
    
    local treemix_dir="${INSTALL_DIR}/treemix"
    mkdir -p "$treemix_dir"
    
    # TreeMix can be installed via conda or compiled
    if command -v conda &> /dev/null; then
        echo "  Installing TreeMix via conda..."
        conda install -c bioconda treemix -y 2>/dev/null && \
            echo -e "${GREEN}  ✓ TreeMix installed via conda${NC}"
    else
        echo -e "${YELLOW}  TreeMix installation requires conda or manual compilation${NC}"
        echo "  See: https://bitbucket.org/nygcresearch/treemix/wiki/Home"
    fi
}

# Update PATH
update_path() {
    echo ""
    echo -e "${BLUE}Updating PATH...${NC}"
    
    local path_additions=""
    
    # Add tool directories to PATH
    [ -d "${INSTALL_DIR}/plink" ] && path_additions="${path_additions}:${INSTALL_DIR}/plink"
    [ -d "${INSTALL_DIR}/admixture" ] && path_additions="${path_additions}:${INSTALL_DIR}/admixture"
    [ -d "${INSTALL_DIR}/eigensoft" ] && path_additions="${path_additions}:${INSTALL_DIR}/eigensoft"
    [ -d "${INSTALL_DIR}/admixtools/bin" ] && path_additions="${path_additions}:${INSTALL_DIR}/admixtools/bin"
    [ -d "${INSTALL_DIR}/rfmix" ] && path_additions="${path_additions}:${INSTALL_DIR}/rfmix"
    [ -d "${INSTALL_DIR}/treemix" ] && path_additions="${path_additions}:${INSTALL_DIR}/treemix"
    
    # Remove leading colon
    path_additions="${path_additions#:}"
    
    if [ -n "$path_additions" ]; then
        echo ""
        echo "Add the following to your ~/.bashrc or ~/.bash_profile:"
        echo ""
        echo -e "${GREEN}export PATH=\"\$PATH:${path_additions}\"${NC}"
        echo ""
        
        # Optionally add to .bashrc
        echo "Would you like to automatically add this to ~/.bashrc? [y/N]"
        read -r response
        if [[ "$response" =~ ^[Yy]$ ]]; then
            echo "" >> ~/.bashrc
            echo "# Ancestry analysis pipeline tools" >> ~/.bashrc
            echo "export PATH=\"\$PATH:${path_additions}\"" >> ~/.bashrc
            echo -e "${GREEN}Added to ~/.bashrc. Run 'source ~/.bashrc' to apply.${NC}"
        fi
    fi
}

# Verify installations
verify_installations() {
    echo ""
    echo -e "${BLUE}Verifying installations...${NC}"
    echo ""
    
    local all_ok=true
    
    # Check PLINK
    if $INSTALL_PLINK; then
        if [ -x "${INSTALL_DIR}/plink/plink" ]; then
            local version
            version=$("${INSTALL_DIR}/plink/plink" --version 2>&1 | head -n1)
            echo -e "${GREEN}  ✓ PLINK: ${version}${NC}"
        else
            echo -e "${RED}  ✗ PLINK: Not found${NC}"
            all_ok=false
        fi
    fi
    
    # Check ADMIXTURE
    if $INSTALL_ADMIXTURE; then
        if [ -x "${INSTALL_DIR}/admixture/admixture" ]; then
            echo -e "${GREEN}  ✓ ADMIXTURE: Installed${NC}"
        else
            echo -e "${RED}  ✗ ADMIXTURE: Not found${NC}"
            all_ok=false
        fi
    fi
    
    # Check EIGENSOFT
    if $INSTALL_EIGENSOFT; then
        if [ -x "${INSTALL_DIR}/eigensoft/smartpca" ]; then
            echo -e "${GREEN}  ✓ EIGENSOFT (smartpca): Installed${NC}"
        else
            echo -e "${YELLOW}  ⚠ EIGENSOFT: Not installed (optional)${NC}"
        fi
    fi
    
    # Check ADMIXTOOLS
    if $INSTALL_ADMIXTOOLS; then
        if [ -x "${INSTALL_DIR}/admixtools/bin/qpAdm" ]; then
            echo -e "${GREEN}  ✓ ADMIXTOOLS (qpAdm): Installed${NC}"
        else
            echo -e "${YELLOW}  ⚠ ADMIXTOOLS: Not installed (optional)${NC}"
        fi
    fi
    
    echo ""
    if $all_ok; then
        echo -e "${GREEN}All core tools installed successfully!${NC}"
    else
        echo -e "${YELLOW}Some tools may need manual installation.${NC}"
    fi
}

# Main
main() {
    print_banner
    detect_os
    
    # Create installation directory
    echo ""
    echo -e "${BLUE}Installation directory: ${INSTALL_DIR}${NC}"
    mkdir -p "$INSTALL_DIR"
    
    # Check prerequisites
    if ! check_prerequisites; then
        echo ""
        echo -e "${YELLOW}Please install missing prerequisites and run this script again.${NC}"
        exit 1
    fi
    
    # Install tools
    if $INSTALL_PLINK; then
        install_plink
    fi
    
    if $INSTALL_ADMIXTURE; then
        install_admixture
    fi
    
    if $INSTALL_EIGENSOFT; then
        install_eigensoft
    fi
    
    if $INSTALL_ADMIXTOOLS; then
        install_admixtools
    fi
    
    if $INSTALL_RFMIX; then
        install_rfmix
    fi
    
    if $INSTALL_TREEMIX; then
        install_treemix
    fi
    
    # Update PATH
    update_path
    
    # Verify
    verify_installations
    
    echo ""
    echo "========================================================"
    echo "                  Installation Complete                  "
    echo "========================================================"
    echo ""
}

# Run main
main
