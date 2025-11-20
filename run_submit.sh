#!/bin/bash
# This script sets up the environment for EIC simulations using EPIC and EICrecon.
# It installs the necessary software, compiles custom configurations, and generates job submission scripts.
# It also handles the simulation and analysis of events based on user-defined parameters.

# example of usage in BASH(!) shell  - not native .csh shell:
# TILE_SIZE=15 N_LAYERS=15 SCINTILLATOR_THICKNESS=0.5 ABSORBER_THICKNESS=5 MOMENTUM=0.5 ./run_submit.sh

# in tcsh shell:
# setenv TILE_SIZE 15
#./run.sh

# The workflow includes:
# 1. Setting up the eic-shell
# 2. Compiling EPIC+EICrecon in your directory if not already done.
# 3. Compiling EPIC+EICrecon with custom specification if different
# 4. Generating a main simulation script = job with ddsim + eicrecon + run root macro.
# 5. Generating a job submission script for HTCondor and submitting it.

#================================================================
#   CONFIGURATION
#================================================================

# detector parameters
DEFAULT_TILE_SIZE=10               # cm
DEFAULT_ABSORBER_THICKNESS=4       # cm
DEFAULT_SCINTILLATOR_THICKNESS=2.4 # cm
DEFAULT_N_LAYERS=10                # 10 layers
DEFAULT_LAYER_GAP=0.1              # cm
MAXIMUM_NHCAL_LENGTH=70            # cm
#calculate default length
DEFAULT_NHCAL_LENGTH=$(echo "$DEFAULT_N_LAYERS * ($DEFAULT_ABSORBER_THICKNESS + $DEFAULT_SCINTILLATOR_THICKNESS + $DEFAULT_LAYER_GAP)" | bc)

# simulation parameters
DEFAULT_MOMENTUM=1          # 1 GeV
DEFAULT_PHI=45              # 45 degrees
DEFAULT_THETA=170           # 170 degrees
DEFAULT_PARTICLE=neutron    # neutron , proton, mu-, pi+
DEFAULT_NUMBER_OF_EVENTS=10 # 10 events per job
DEFAULT_JOBS=10             # 10 jobs

# user-defined values - can be set via command line or environment variables
TILE_SIZE=${TILE_SIZE:-$DEFAULT_TILE_SIZE}
ABSORBER_THICKNESS=${ABSORBER_THICKNESS:-$DEFAULT_ABSORBER_THICKNESS}
SCINTILLATOR_THICKNESS=${SCINTILLATOR_THICKNESS:-$DEFAULT_SCINTILLATOR_THICKNESS}
LAYER_GAP=${LAYER_GAP:-$DEFAULT_LAYER_GAP}
N_LAYERS=${N_LAYERS:-$DEFAULT_N_LAYERS}

# simulation parameters - can be set via command line or environment variables
MOMENTUM=${MOMENTUM:-$DEFAULT_MOMENTUM}
PHI=${PHI:-$DEFAULT_PHI}
THETA=${THETA:-$DEFAULT_THETA}
PARTICLE=${PARTICLE:-$DEFAULT_PARTICLE}
NUMBER_OF_EVENTS=${NUMBER_OF_EVENTS:-$DEFAULT_NUMBER_OF_EVENTS}
JOBS=${JOBS:-$DEFAULT_JOBS}

#================================================================
DETECTOR_CONFIG="epic_backward_hcal_only.xml" # "epic_full.xml" or "epic_backward_hcal_only.xml"
#  make unique dir name for whole configuration
NHCAL_CONFIG="nhcal_only_tile${TILE_SIZE}cm_absorber${ABSORBER_THICKNESS}cm_scintillator${SCINTILLATOR_THICKNESS}cm_${N_LAYERS}layers"
SIM_CONFIG="${PARTICLE}_p${MOMENTUM}gev_phi${PHI}_theta${THETA}_${NUMBER_OF_EVENTS}events"

# output directories
output_dir="/gpfs02/eic/${USER}/output/${NHCAL_CONFIG}/${SIM_CONFIG}"
my_epic_dir="/gpfs02/eic/${USER}/epic"
my_eicrecon_dir="/gpfs02/eic/${USER}/EICrecon"
current_dir="$(pwd)"

#================================================================
#   UTILITY FUNCTIONS
#================================================================
set -euo pipefail # Exit on error, undefined variables, and pipe failures
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}
error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
    echo "  Context: Function ${FUNCNAME[1]}, Line ${BASH_LINENO[0]}" >&2
}
display_configuration() {
    cat <<EOF

===  SIMULATION CONFIGURATION ===
Detector Parameters:
  - Tile Size: ${TILE_SIZE} cm
  - Absorber Thickness: ${ABSORBER_THICKNESS} cm  
  - Scintillator Thickness: ${SCINTILLATOR_THICKNESS} cm
  - Number of Layers: ${N_LAYERS}
  - Calculated Length: ${NEW_LENGTH} cm

Simulation Parameters:
  - Particle: ${PARTICLE}
  - Momentum: ${MOMENTUM} GeV
  - Phi: ${PHI}°, Theta: ${THETA}°
  - Events per Job: ${NUMBER_OF_EVENTS}
  - Number of Jobs: ${JOBS}

Output Directory: ${output_dir}
===================================
EOF
}
#=================================================================
# 1. Setting up the eic-shell
# Setup EIC shell environment
# Globals:
#   EICSHELL - Path to eic-shell executable
#   USER - Current username
# Returns:
#   0 on success, 1 on failure
#================================================================
setup_eicshell() {
    export EICSHELL="/eic/u/${USER}/eic/eic-shell"
    if [[ ! -f "$EICSHELL" ]]; then
        log "Installing eic-shell..."
        mkdir -p "$HOME/eic"
        cd "$HOME/eic"
        if ! curl -L https://github.com/eic/eic-shell/raw/main/install.sh | bash; then
            error "Failed to install eic-shell"
            exit 1
        fi
        cd "$current_dir"
    fi
}
#=================================================================

#=================================================================
# 2. Compiling EPIC in your directory if not already done.
#=================================================================

setup_epic() {
    log "Checking EPIC setup..."
    # Check if EPIC is already properly installed
    if [[ -d "$my_epic_dir" ]] && [[ -f "$my_epic_dir/install/bin/thisepic.sh" ]]; then
        log "EPIC already installed at $my_epic_dir"
        return 0
    fi

    log "Setting up EPIC..."

    # Create parent directory if it doesn't exist
    local parent_dir="$(dirname "$my_epic_dir")"
    mkdir -p "$parent_dir"

    # Remove any incomplete installation
    if [[ -d "$my_epic_dir" ]]; then
        log "Removing incomplete EPIC installation..."
        rm -rf "$my_epic_dir"
    fi

    # Clone EPIC
    cd "$parent_dir"
    log "Cloning EPIC repository..."
    if ! git clone https://github.com/eic/epic.git "$(basename "$my_epic_dir")"; then
        error "Failed to clone EPIC repository"
        exit 1
    fi

    cd "$my_epic_dir"
    # Build EPIC

    cat <<EOF | $EICSHELL
    echo "Configuring EPIC build..."
    # Use ccache for faster builds
    export CC="ccache gcc"
    export CXX="ccache g++"

    # Run CMake configure step with ccache launcher support
    echo " Building and installing EPIC (this may take several minutes)..."
    if ! cmake -B build -S . \
        -DCMAKE_C_COMPILER_LAUNCHER=ccache \
        -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
        -DCMAKE_INSTALL_PREFIX=install; then
        echo "CMake configuration failed."
        exit 1
    fi
    if ! cmake --build build -j$(nproc) -- install; then
        exit 1
    fi
EOF
    # Verify installation
    if [[ ! -f "$my_epic_dir/install/bin/thisepic.sh" ]]; then
        error "EPIC installation incomplete - missing thisepic.sh"
        exit 1
    fi

    cd "$current_dir"
    log "EPIC setup completed successfully"

}

#================================================================
#   3. Compiling EICrecon in your directory if not already done.
#================================================================

setup_eicrecon() {
    log "Checking EICrecon setup..."
    # First, ensure EPIC is properly installed
    if [[ ! -f "$my_epic_dir/install/bin/thisepic.sh" ]]; then
        log "Expected EPIC installation at: $my_epic_dir/install/bin/thisepic.sh"
        error "EPIC installation not found. Please set up EPIC first."

        exit 1
    fi
    # Check if EICrecon is already properly installed
    if [[ -d "$my_eicrecon_dir" ]] && [[ -f "$my_eicrecon_dir/install/bin/eicrecon-this.sh" ]]; then
        log "EICrecon already installed at $my_eicrecon_dir"
        return 0
    fi

    log "Setting up EICrecon..."

    # Create parent directory if it doesn't exist
    local parent_dir="$(dirname "$my_eicrecon_dir")"
    mkdir -p "$parent_dir"

    # Remove any incomplete installation
    if [[ -d "$my_eicrecon_dir" ]]; then
        log "Removing incomplete EICrecon installation..."
        rm -rf "$my_eicrecon_dir"
    fi

    # Clone EICrecon
    cd "$parent_dir"
    log "Cloning EICrecon repository..."
    if ! git clone https://github.com/eic/EICrecon.git "$(basename "$my_eicrecon_dir")"; then
        error "Failed to clone EICrecon repository"
        exit 1
    fi

    # Build EICrecon
    cd "$my_eicrecon_dir"

    log "Configuring EICrecon build..."
    cat <<EOF | $EICSHELL
    source "$my_epic_dir/install/bin/thisepic.sh"
    if ! cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install; then
        exit 1
    fi
    echo "Building and installing EICrecon (this may take several minutes)..."
    if ! cmake --build build -j$(nproc) -- install; then
        exit 1
    fi
    echo "Installing EICrecon..."
    if ! cmake --install build; then
        exit 1
    fi
EOF
    # Verify installation
    if [[ ! -f "$my_eicrecon_dir/install/bin/eicrecon-this.sh" ]]; then
        error "EICrecon installation incomplete - missing eicrecon-this.sh"
        exit 1
    fi

    cd "$current_dir"
    log "EICrecon setup completed successfully"
}
#=================================================================

#================================================================
#   3. Compiling EPIC+EICrecon with custom specification if different
#================================================================

compile_custom_epic() {

    # check if any detector parameters are different from default
    if [[ "$TILE_SIZE" == "$DEFAULT_TILE_SIZE" ]] &&
        [[ "$ABSORBER_THICKNESS" == "$DEFAULT_ABSORBER_THICKNESS" ]] &&
        [[ "$SCINTILLATOR_THICKNESS" == "$DEFAULT_SCINTILLATOR_THICKNESS" ]] &&
        [[ "$N_LAYERS" == "$DEFAULT_N_LAYERS" ]]; then
        log "No custom EPIC compilation needed - using default parameters"
        return 0
    fi

    # Create working copies
    local work_epic_dir="$output_dir/../epic"
    local work_eicrecon_dir="$output_dir/../EICrecon"

    # check if EPIC has already been compiled with custom parameters in the output directory
    if [[ -d "$work_epic_dir" ]] && [[ -f "$work_epic_dir/install/bin/thisepic.sh" ]]; then
        log "Custom EPIC already compiled at $work_epic_dir"
        return 0
    fi

    log "Compiling EPIC with custom specifications..."
    log "Tile size: $TILE_SIZE cm, Absorber thickness: $ABSORBER_THICKNESS cm, Scintillator thickness: $SCINTILLATOR_THICKNESS cm, Number of layers: $N_LAYERS"

    log "Copying EPIC to working directory..."

    # copy without build directory, .git and install
    rsync -av --exclude build/ --exclude install/ --exclude .git/ "$my_epic_dir/" "$work_epic_dir/"

    log "Modifying EPIC configuration..."
    # Modify EPIC configuration
    local nhcal_config="$work_epic_dir/compact/hcal/backward_template.xml"
    local epic_definitions="$work_epic_dir/compact/definitions.xml"

    # More robust sed patterns with backup
    cp "$nhcal_config" "${nhcal_config}.backup"
    cp "$epic_definitions" "${epic_definitions}.backup"

    # check if TILE_SIZE is different from default

    if [[ "$TILE_SIZE" != "$DEFAULT_TILE_SIZE" ]]; then
        local new_size="${TILE_SIZE}0 * mm"
        if sed -i "s/grid_size_x=\"[^\"]*\"/grid_size_x=\"$new_size\"/" "$nhcal_config" &&
            sed -i "s/grid_size_y=\"[^\"]*\"/grid_size_y=\"$new_size\"/" "$nhcal_config"; then #
            log "Tile size updated to $new_size"
        else
            error "Failed to update tile size in $nhcal_config"
            exit 1
        fi
    fi
    #  with the new number of layers it is tricky since it depends on thickness and total length
    # HcalEndcapNSingleLayerThickness=HcalEndcapNSteelThickness + HcalEndcapNPolystyreneThickness + HcalEndcapNLayerGap
    # HcalEndcapNLayer_NRepeat = floor(HcalEndcapN_length / HcalEndcapNSingleLayerThickness)
    # check if ABSORBER_THICKNESS is different from default
    if [[ "$ABSORBER_THICKNESS" != "$DEFAULT_ABSORBER_THICKNESS" ]]; then
        if sed -i "s/\(HcalEndcapNSteelThickness\"[[:space:]]*value=\"\)[^ ]*/\1$ABSORBER_THICKNESS/" "$nhcal_config"; then
            log "Absorber thickness updated to $ABSORBER_THICKNESS"
        else
            error "Failed to update absorber thickness in $nhcal_config"
            exit 1
        fi
    fi

    # check if SCINTILLATOR_THICKNESS is different from default
    if [[ "$SCINTILLATOR_THICKNESS" != "$DEFAULT_SCINTILLATOR_THICKNESS" ]]; then
        if sed -i "s/\(HcalEndcapNPolystyreneThickness\"[[:space:]]*value=\"\)[^ ]*/\1$SCINTILLATOR_THICKNESS/" "$nhcal_config"; then
            log "Scintillator thickness updated to $SCINTILLATOR_THICKNESS"
        else
            error "Failed to update scintillator thickness in $nhcal_config"
            exit 1
        fi
    fi
    # change the total length of nHcal to the new length if it is different from default
    if [[ "$NEW_LENGTH" != "$DEFAULT_NHCAL_LENGTH" ]]; then
        if sed -i "s/\(HcalEndcapN_length\"[[:space:]]*value=\"\)[^*]*/\1$NEW_LENGTH/" "$epic_definitions"; then
            log "NHCAL length updated to $NEW_LENGTH"
        else
            error "Failed to update NHCAL length in $epic_definitions"
            exit 1
        fi
    fi

    log "Compiling custom EPIC..."

    cat <<EOF | $EICSHELL
    cd "$work_epic_dir"

    export CC="ccache gcc"
    export CXX="ccache g++"

    if ! cmake -B build -S . \
        -DCMAKE_C_COMPILER_LAUNCHER=ccache \
        -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
        -DCMAKE_INSTALL_PREFIX=install; then
        echo "Failed to configure EPIC build" >&2
        exit 1
    fi
    if ! cmake --build build -j8 -- install; then
        echo "Failed to build EPIC" >&2
        exit 1
    fi
EOF
}
#=================================================================

#================================================================
#   MAIN SCRIPT GENERATION
#================================================================

generate_main_script() {
    log "Generating main simulation script..."

    generated_script="generated_script_$NHCAL_CONFIG.sh"

    cat >$generated_script <<'EOF'
#!/bin/bash
cd OUTPUT_DIR_PLACEHOLDER
# Parse arguments
if [[ $# -ne 7 ]]; then
    echo "Usage: $0 <Cluster> <Process> <Momentum> <Phi> <Theta> <Particle> <NumberOfEvents>"
    echo "current values: $CLUSTER $PROCESS $GUN_MOMENTUM $GUN_PHI $GUN_THETA $PARTICLE $NUMBER_OF_EVENTS"
    exit 1
fi
export CLUSTER="$1"
export PROCESS="$2"
export GUN_MOMENTUM="$3"
export GUN_PHI="$4"
export GUN_THETA="$5"
export PARTICLE="$6"
export NUMBER_OF_EVENTS="$7"

# Set up file names
export FILENAME="${PARTICLE}_${NUMBER_OF_EVENTS}events_p${GUN_MOMENTUM}gev_phi${GUN_PHI}_theta${GUN_THETA}_job${CLUSTER}_${PROCESS}"
export DDSIM_FILE="sim_${FILENAME}.edm4hep.root"
export EICRECON_FILE="eicrecon_${FILENAME}.edm4eic.root"

# Calculate angle ranges (small ranges for single-angle shooting)
export GUN_THETA_MIN=$(echo "$GUN_THETA - 0.0001" | bc -l)
export GUN_THETA_MAX=$(echo "$GUN_THETA + 0.0001" | bc -l)
export GUN_PHI_MIN=$(echo "$GUN_PHI - 0.0001" | bc -l)
export GUN_PHI_MAX=$(echo "$GUN_PHI + 0.0001" | bc -l)
export GUN_MOMENTUM_MIN=$(echo "$GUN_MOMENTUM - 0.00001" | bc -l)
export GUN_MOMENTUM_MAX=$(echo "$GUN_MOMENTUM + 0.00001" | bc -l)

cat << 'EOFINNER' | EICSHELL_PLACEHOLDER
    # Source environment
    if [[ -f "OUTPUT_DIR_PLACEHOLDER/epic/install/bin/thisepic.sh" ]]; then
        source "OUTPUT_DIR_PLACEHOLDER/epic/install/bin/thisepic.sh" epic
    elif [[ -f "MY_EPIC_DIR_PLACEHOLDER/install/bin/thisepic.sh" ]]; then
        source "MY_EPIC_DIR_PLACEHOLDER/install/bin/thisepic.sh" epic
    else
        error "EPIC installation not found"
        exit 1
    fi

    if ! ddsim \
            --compactFile "$DETECTOR_PATH/DETECTOR_CONFIG_PLACEHOLDER" \
            --numberOfEvents "$NUMBER_OF_EVENTS" \
            --random.seed "$(date +%N)" \
            --enableGun \
            --gun.particle "$PARTICLE" \
            --gun.thetaMin "${GUN_THETA_MIN}*degree" \
            --gun.thetaMax "${GUN_THETA_MAX}*degree" \
            --gun.phiMin "${GUN_PHI_MIN}*degree" \
            --gun.phiMax "${GUN_PHI_MAX}*degree" \
            --gun.distribution uniform \
            --gun.momentumMin "${GUN_MOMENTUM_MIN}*GeV" \
            --gun.momentumMax "${GUN_MOMENTUM_MAX}*GeV" \
            --outputFile "$DDSIM_FILE"; then
            echo "DDSIM simulation failed"
            exit 1
    fi

    # Source EICrecon
    if [[ -f "MY_EICRECON_DIR_PLACEHOLDER/install/bin/eicrecon-this.sh" ]]; then
        source "MY_EICRECON_DIR_PLACEHOLDER/install/bin/eicrecon-this.sh" epic
    else
        echo "EICrecon installation not found"
        exit 1
    fi

    # Run EICrecon if needed
    echo "Running EICrecon..."
    if ! eicrecon "$DDSIM_FILE" \
            -Ppodio:output_file="$EICRECON_FILE"
            -Ppodio:output_collections="MCParticles,HcalEndcapNClusters,HcalEndcapNTruthClusters,\
        EcalEndcapNClusters,EcalEndcapNTruthClusters,\
        HcalBarrelClusters,HcalBarrelTruthClusters"; then
        echo "EICrecon failed"
        exit 1
    fi

    # Run analysis
    echo "Running ROOT analysis..."
    analysis_script="CURRENT_DIR_PLACEHOLDER/example_macro.C"
    if [[ ! -f "$analysis_script" ]]; then
        echo "Analysis script not found: $analysis_script"
        exit 1
    fi
    output_file="ana_${FILENAME}.root"
    if ! root -l -b -q "${analysis_script}(\\\"${EICRECON_FILE}\\\", \\\"${output_file}\\\")"; then
        echo "ROOT analysis failed"
        exit 1
    fi
    echo "Job completed successfully"
EOFINNER
EOF
    # Replace placeholders with actual values
    sed -i "s|OUTPUT_DIR_PLACEHOLDER|$output_dir|g" $generated_script
    sed -i "s|MY_EPIC_DIR_PLACEHOLDER|$my_epic_dir|g" $generated_script
    sed -i "s|MY_EICRECON_DIR_PLACEHOLDER|$my_eicrecon_dir|g" $generated_script
    sed -i "s|CURRENT_DIR_PLACEHOLDER|$current_dir|g" $generated_script
    sed -i "s|EICSHELL_PLACEHOLDER|$EICSHELL|g" $generated_script
    sed -i "s|DETECTOR_CONFIG_PLACEHOLDER|$DETECTOR_CONFIG|g" $generated_script

    chmod +x $generated_script

    log "Main simulation script generated successfully"
}

#================================================================
#   JOB SUBMISSION SCRIPT
#================================================================

generate_job_script() {
    local temp_job="$current_dir/generated.job"
    # Remove existing job file if it exists
    if [[ -f "$temp_job" ]]; then
        rm -f "$temp_job"
    fi
    cat >"$temp_job" <<EOF
Universe                = vanilla
GetEnv                  = False
Requirements            = (CPU_Speed >= 1)
Rank                    = CPU_Speed
Initialdir              = $current_dir
Arguments               = \$(Cluster) \$(Process) $MOMENTUM $PHI $THETA $PARTICLE $NUMBER_OF_EVENTS
Executable              = $generated_script
Error                   = $output_dir/log/error\$(Cluster)_\$(Process).err
Output                  = $output_dir/log/out\$(Cluster)_\$(Process).out
Log                     = $output_dir/log/log\$(Cluster)_\$(Process).log
# File transfer settings
Should_Transfer_Files   = YES
When_To_Transfer_Output = ON_EXIT
Transfer_Input_Files    = $generated_script
Transfer_Output_Files   = ""

Queue $JOBS
EOF
    if [[ ! -f "$temp_job" ]]; then
        return 1
    fi
    # Only echo the filename - no log messages
    echo "$temp_job"
}

#================================================================
#   MAIN EXECUTION
#================================================================

main() {
    log "=== Starting EPIC Simulation ==="
    #================================================================
    # first is to calculate new length
    NEW_LENGTH=$(echo "$N_LAYERS * ($ABSORBER_THICKNESS + $SCINTILLATOR_THICKNESS + $LAYER_GAP)" | bc)
    # check if new length is larger than maximum
    if (($(echo "$NEW_LENGTH > $MAXIMUM_NHCAL_LENGTH" | bc -l))); then
        error "New length $NEW_LENGTH cm exceeds maximum length of $MAXIMUM_NHCAL_LENGTH cm"
        exit 1
    fi
    #================================================================
    display_configuration

    # Setup phase
    log "Setting up eic-shell..."
    setup_eicshell
    setup_epic
    setup_eicrecon

    # Create output directories
    mkdir -p "$output_dir/log"
    rm -f $output_dir/log/*.*
    log "Output directory created: $output_dir"

    compile_custom_epic

    generate_main_script

    log "Generating HTCondor job submission script..."
    local job_file=$(generate_job_script)

    if [[ -z "$job_file" ]] || [[ ! -f "$job_file" ]]; then
        error "Failed to create job submission file"
        exit 1
    fi
    log "Job submission file created: $job_file"
    log "Setup completed successfully!"

    log "Submitting job..."
    if condor_submit "$job_file"; then
        log "Job submitted successfully! Output files will be saved in $output_dir"
    else
        error "Job submission failed"
        exit 1
    fi
    # wait for the jobs to finish
    #log "Waiting for jobs to finish..."
    # ./condor_control.sh
    log "=== EPIC Simulation Completed ==="
    #  now one can merge the output files into one
    # hadd -f -j -k "${output_dir}/merged_ana.root" "${output_dir}/ana*.root"
}
# Run main function
main "$@"
