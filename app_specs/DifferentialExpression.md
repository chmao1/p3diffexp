
# Application specification: DifferentialExpression

This is the application specification for service with identifier DifferentialExpression.

The backend script implementing the application is [App-DifferentialExpression.pl](../service-scripts/App-DifferentialExpression.pl).

The raw JSON file for this specification is [DifferentialExpression.json](DifferentialExpression.json).

This service performs the following task:   Parses and transforms users differential expression data

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| xfile | Experiment Data File | WS: ExpList  | :heavy_check_mark: |  |
| mfile | Metadata File | WS: ExpMetadata  |  |  |
| ustring | User string | string  | :heavy_check_mark: |  |
| output_path | Output Folder | folder  |  |  |
| output_file | File Basename | wsid  |  |  |

