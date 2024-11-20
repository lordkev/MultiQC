from typing import Dict, Union
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table
import logging

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Yleaf",
            anchor="yleaf",
            href="https://github.com/genid/Yleaf",
            info="Yleaf: software for human Y-chromosomal haplogroup inference from next generation sequencing data",
            doi=["10.1093/molbev/msy032"],
        )
        yleaf_data: Dict = dict()

        for f in self.find_log_files("yleaf"):
            parsed_data = self.parse_logs(f)
            yleaf_data.update(parsed_data)

        yleaf_data = self.ignore_samples(yleaf_data)

        if len(yleaf_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(yleaf)} Yleaf reports")

        # Write data to file
        self.write_data_file(yleaf, "multiqc_yleaf")

        self.add_software_version(None)

        # Headers dictionary for Contamination Data
        headers = {
            "Hg": {
                "title": "Haplogroup",
                "description": "Y-chromosomal haplogroup.",
                "scale": False,
            },
            "Hg_marker": {
                "title": "Haplogroup Marker",
                "description": "Y-chromosomal haplogroup marker.",
                "scale": False,
            },
            "Hg_marker": {
                "title": "Haplogroup Marker",
                "description": "Y-chromosomal haplogroup marker.",
                "scale": False,
            },
            "Total_reads": {
                "title": "Total Reads",
                "description": "Total number of reads.",
                "min": 0,
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "Valid_markers": {
                "title": "Valid Markers",
                "description": "Number of valid markers.",
                "min": 0,
                "scale": "Blues",
                "format": "{:,.0f}",
            },
            "QC-score": {
                "title": "QC Score",
                "description": "Quality control score.",
                "min": 0,
                "scale": "Oranges",
                "format": "{:,.2f}",
            },
            "QC-1": {
                "title": "QC-1",
                "description": "Quality control 1.",
                "min": 0,
                "scale": "Oranges",
                "format": "{:,.2f}",
            },
            "QC-2": {
                "title": "QC-2",
                "description": "Quality control 2.",
                "min": 0,
                "scale": "Oranges",
                "format": "{:,.2f}",
            },
            "QC-3": {
                "title": "QC-3",
                "description": "Quality control 3.",
                "min": 0,
                "scale": "Oranges",
                "format": "{:,.2f}",
            },
        }

        self.add_section(
            name="Yleaf",
            anchor="yleaf-section",
            description="Yleaf: software for human Y-chromosomal haplogroup inference from next generation sequencing data",
            plot=table.plot(
                data=yleaf_data,
                headers=headers,
                pconfig={
                    "id": "yleaf-table",
                    "title": "Yleaf",
                },
            ),
        )

        for header in headers.values():
            header["hidden"] = True
        headers["Contamination Status"]["hidden"] = False
        headers["Hg"]["hidden"] = False
        headers["Valid_markers"]["hidden"] = False
        headers["QC-score"]["hidden"] = False

        self.general_stats_addcols(yleaf_data, headers)

    def parse_logs(self, f: str) -> Dict[str, Dict[str, Union[float, str]]]:
        parsed_data: Dict[str, Dict[str, Union[float, str]]] = {}
        file_content = f["f"]
        lines = file_content.strip().splitlines()

        # Assuming the first line contains headers
        headers = [header.strip().replace('"', "") for header in lines[0].strip().split("\t")[1:]]

        for line in lines[1:]:
            columns = line.strip().split("\t")
            sample_name = columns[0]  # First column is the sample name
            s_name = sample_name.replace('"', "")  # Remove quotes
            values = columns[1:]
            self.add_data_source(f, s_name)

            # Map headers to values for this sample
            parsed_data[s_name] = {
                key: (
                    float(value.strip().replace('"', ""))
                    if value.strip().replace('"', "").replace(".", "", 1).isdigit()
                    else value.strip().replace('"', "")
                )
                for key, value in zip(headers, values)
            }

        return parsed_data
