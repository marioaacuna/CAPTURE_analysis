# Extended CAPTURE Analysis Pipeline

This extended version of the CAPTURE analysis pipeline is designed to handle multi-condition experiments with miniscope recordings. It maintains the core functionality of the original pipeline while adding support for more complex experimental designs.

## Key Features

- Support for multiple experimental conditions:
  - Baseline recordings
  - Saline/Formalin injections (paired design)
  - Sham/SNI surgeries
- Integration of both behavioral (6-cam) and neural (miniscope) data
- Structured metadata management
- Maintained compatibility with original analysis methods

## Directory Structure

```
extended_pipeline/
├── configs/
│   └── metadata.json       # Experiment metadata configuration
├── Analysis_MA/           # Motion analysis scripts
├── Behavioral_analysis/   # Behavioral data processing
├── Preprocessing/        # Data preprocessing scripts
├── Utility/             # Helper functions
└── VideoAnalysis/       # Video processing tools
```

## Metadata Configuration

The pipeline uses a JSON configuration file (`configs/metadata.json`) to manage experimental metadata. This file contains:
- Data root directory paths
- Animal IDs and their conditions
- Recording types (6-cam/miniscope availability)

Example metadata structure:
```json
{
    "data_root_dir": "path/to/data",
    "conditions": {
        "baseline": {
            "data_path": "baseline",
            "animals": [
                {
                    "id": "1234",
                    "has_6cam": true,
                    "has_miniscope": true,
                    "path": "ID_1234"
                }
            ]
        }
    }
}
```

## Data Organization

The pipeline expects data to be organized as follows:
```
root_directory/
├── baseline/
│   ├── 6cam_data/
│   │   └── ID_XXXX/
│   │       └── YYYYMMDD/
│   │           └── DANNCE/
│   └── miniscope_data/
├── Formalin_injection/
├── saline_injection/
├── sham/
└── SNI/
```

## Key Scripts

1. `general_configs.m`: Configuration settings for the extended pipeline
2. Script modifications for condition handling:
   - Automatic date folder detection
   - Multi-condition data loading
   - Integrated metadata parsing

## Usage

1. Configure metadata:
   - Update `metadata.json` with your experimental structure
   - Verify all paths and animal IDs

2. Run analysis:
   ```matlab
   % Basic initialization
   close all force
   clc
   clear
   cd('path/to/CAPTURE_analysis')
   
   % Add necessary paths
   addpath('3rd party toolboxes')
   addpath('Common')
   addpath('MotionMapper-master')
   addpath('Species_specific_files')
   addpath('Clustering')
   addpath('Animating')
   addpath(genpath('extended_pipeline'))
   ```

## Key Differences from Original Pipeline

1. **Metadata Management**:
   - JSON-based configuration instead of hardcoded paths
   - Structured experimental condition handling

2. **Data Organization**:
   - Support for multiple experimental conditions
   - Automated date folder handling
   - Integration of behavioral and neural data

3. **Analysis Flow**:
   - Condition-aware data loading
   - Support for paired experimental designs
   - Enhanced error handling for missing data

## Contributing

When modifying this pipeline:
1. Maintain compatibility with original analysis methods
2. Document changes in metadata structure
3. Update error handling for new data structures
4. Test across all experimental conditions

## Dependencies

Same as original pipeline, plus:
- MATLAB JSON parsing capability (built-in)
- Original CAPTURE analysis toolboxes

## Notes

- Always verify metadata configuration before running analysis
- Check data paths and structure match expected organization
- Backup analysis results for each experimental condition separately