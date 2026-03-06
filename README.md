# LT SWAT Point Source (PS) Data Processing Toolbox

This repository contains **R scripts** and tools for processing point source (PS) discharge data—such as wastewater treatment plant outputs—to prepare inputs for the **SWAT (Soil and Water Assessment Tool)** model, specifically tailored for Lithuania's river modeling system (LRMS). The workflow focuses on data validation, pollutant load calculation, and risk assessment for water bodies.

## Purpose

The purpose of this toolbox is to **automate the integration of point source pollution data into the LT SWAT modeling framework** and to support the **assessment of the impacts of point sources on water bodies**. This is part of the broader **LIFE22-IPE-LT-LIFE-SIP-Vanduo** project, which aims to enhance integrated water management in Lithuania.

## Project Overview

- **Created on**: 2025-06-05  
- **Last Modified**: 2026-03-05  
- **Author**: Svajunas Plunge  
- **Email**: svajunas_plunge@sggw.edu.pl  

### Key Features
1.  **Quality Assurance (QA)**: The `ps_qa.R` script is used for data validation and integrity checks to ensure point source datasets are  correctly represented in LRMS.
2.  **Impact Analysis**: `ps_impacts.R` calculates pollutant loads and evaluates their specific impacts on the water bodies (WB).
3.  **Risk Assessment**: `ps_risk_wb.R` identifies risk WBs comparing *base* and *change* scenarios. 
4.  **Modular Functions**: Core logic and shared utilities are maintained in `function.R` to ensure consistency across different analysis steps.

### Repository Structure
```
LT-SWAT-PS-Processing/
├── Data/                # Directory for storing raw and processed datasets
├── function.R           # Collection of custom R functions used across scripts
├── ps_impacts.R         # Script for assessing environmental impacts of point sources
├── ps_qa.R              # Quality assurance and data cleaning for point source datasets
├── ps_risk_wb.R         # Water body risk assessment based on point source loads
├── update_small.R       # Utility script for updating small point source datasets
├── .Renviron            # Environment variables configuration (e.g., Postgres database connection)
├── ps_in_ltswat.Rproj   # RStudio project file for the LT-SWAT PS processing workflow
└── README.md            # Project documentation and usage instructions
```

## Prerequisites

### Software Requirements
*   **R** (100% of the project is written in R).
*   **RStudio** is recommended for managing the `.Rproj` environment.

### Data Requirements
*   Point source discharge datasets as required by the analysis scripts.
*   A properly configured `.Renviron` file to define local data directories and environment settings.
*   Access to a PostGres database if required for data storage and retrieval.
*   LRMS (Lithuanian River Modeling System) data for integrating point source information.

### Usage
1.  **Clone the Repository**: Download or clone the `ps_in_ltswat` repository from GitHub.
2.  **Initialize Project**: Open `ps_in_ltswat.Rproj` in RStudio to set the working directory.
3.  **Configure Environment**: Edit the `.Renviron` file to specify your local file paths.
4.  **Load Functions**: Source `function.R` to load the necessary utilities for the project.
5.  **Run Scripts**:
    *   Execute `ps_qa.R` to validate your data.
    *   Execute `ps_impacts.R` or `ps_risk_wb.R` to generate model inputs and risk reports.

## Acknowledgments

This work was carried out within the **LIFE22-IPE-LT-LIFE-SIP-Vanduo** project (*Integrated water management in Lithuania*, ref: LIFE22-IPE-LT-LIFE-SIP-Vanduo/101104645), funded by the **European Union LIFE program** under the grant agreement No 101104645.
