i# Setting up RDG with Pixi

This guide will help you set up the RDG (Ribosome Decision Graph) project using pixi for dependency management.

## Prerequisites

1. **Install pixi** if you haven't already:
   ```bash
   # On macOS/Linux
   curl -fsSL https://pixi.sh/install.sh | bash
   
   # On Windows (PowerShell)
   iwr -useb https://pixi.sh/install.ps1 | iex
   ```

2. **Restart your terminal** or source your shell configuration to make pixi available.

## Setup Steps

1. **Initialize the project** (if starting fresh):
   ```bash
   cd /path/to/your/rdg/project
   ```

2. **Install dependencies**:
   ```bash
   pixi install
   ```

3. **Install the package in development mode**:
   ```bash
   pixi run install-dev
   ```

## Available Commands

### Development Commands
```bash
# Run tests
pixi run test

# Run tests with coverage
pixi run test-coverage

# Install package in development mode
pixi run install-dev
```

### CLI Commands
```bash
# Get help for the main CLI
pixi run rdg-help

# Get help for specific subcommands
pixi run rdg-visualise-help
pixi run rdg-construct-help
pixi run rdg-prune-help
```

### Documentation (optional feature)
```bash
# Switch to docs environment and build documentation
pixi shell -e docs
pixi run docs-build

# Clean documentation build
pixi run docs-clean
```

## Usage Examples

Once installed, you can use the RDG CLI directly:

```bash
# Activate the pixi environment
pixi shell

# Now you can use the rdg command directly
rdg --help

# Visualize graphs from a FASTA file
rdg visualise input.fasta --num_visualised 3

# Construct graphs from different file formats
rdg construct input.fasta --input_format fasta
rdg construct input.gtf --input_format gtf
rdg construct input.bed --input_format bed

# Prune graphs based on riboseq data
rdg prune graph.sqlite riboseq.bedgraph --min_read_support 10
```

## Project Structure

The project follows this structure:
```
RDG/
├── pixi.toml              # Pixi configuration
├── requirements.txt       # Python dependencies
├── setup.py              # Package setup
├── RDG/                  # Main package
│   ├── __init__.py
│   ├── cli.py            # Command-line interface
│   ├── RDG.py           # Core RDG classes
│   ├── plot.py          # Plotting functionality
│   ├── RDG_to_file.py   # File I/O operations
│   └── sequence_to_RDG.py # Sequence processing
├── tests/               # Test suite
└── docs/               # Documentation
```

## Features

- **Dependency Management**: All Python and system dependencies managed by pixi
- **Cross-platform**: Works on Linux, macOS, and Windows
- **Reproducible**: Locked dependency versions ensure consistent environments
- **Task Management**: Pre-defined tasks for common operations
- **Multiple Environments**: Separate environments for development and documentation

## Troubleshooting

1. **If dependencies fail to install**: Check that you have access to the conda-forge and bioconda channels
2. **If the CLI command is not found**: Make sure you've run `pixi run install-dev` and are in the pixi shell (`pixi shell`)
3. **Import errors**: Ensure the package is installed in development mode with `pixi run install-dev`

## Next Steps

1. **Add test data**: Place test FASTA, GTF, or BED files in a `test_data/` directory
2. **Customize tasks**: Modify the `[tasks]` section in `pixi.toml` to add project-specific commands
3. **Add CI/CD**: The pixi configuration can be used in GitHub Actions or other CI systems
