# Licensed with the 3-clause BSD license.  See LICENSE for details.

__all__ = ["CSSChecker"]

import os
import re
import email
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Union
from xml.etree.ElementTree import ElementTree

import requests
from pds4_tools import pds4_read
from pds4_tools.reader.label_objects import Label as PDS4Label
from astropy.time import Time

from sbsearch import SBSearch
from sbsearch.logging import ProgressTriangle
from .model import (Observation, Found, Ephemeris, SyncStatus,
                    CatalinaBigelow, CatalinaKittPeak, CatalinaLemmon)


class CSSChecker(SBSearch):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, logger_name='CSS-Checker', **kwargs)

    # URL of most recent file list
    LATEST_FILES: str = (
        "https://sbnarchive.psi.edu/pds4/surveys/catalina_extras/file_list.latest.txt"
    )

    # URL prefix for the CSS archive at PSI
    ARCHIVE_PREFIX: str = "https://sbnarchive.psi.edu/pds4/surveys/"

    # local file list
    LIST_FILE: str = "css-file-list.txt"

    def download_list_file(self) -> None:
        """Check PDS archive for new data and save to file."""

        download: bool = False

        if os.path.exists(self.LIST_FILE):
            # file exists, check for an update
            last_download_date: datetime = datetime.fromtimestamp(
                os.stat(self.LIST_FILE).st_mtime)
            response: requests.Response = requests.head(self.LATEST_FILES)
            try:
                file_date: datetime = datetime(
                    *email.utils.parsedate(response.headers['Last-Modified'])
                )
                if last_download_date < file_date:
                    download = True
                    self.logger.info('New file list available.')
            except KeyError:
                pass
        else:
            # file does not exist, download new file
            download = True

        if download:
            with requests.get(self.LATEST_FILES, stream=True) as r:
                r.raise_for_status()
                with open(self.LIST_FILE, "wb") as f:
                    chunk: bytearray
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
                self.logger.info('Downloaded file list.')

                stat: os.stat_result = os.stat(self.LIST_FILE)
                file_date = Time(stat.st_mtime, format='unix')
                self.logger.info(f"  Size: {stat.st_size / 1048576:.2f} MiB")
                self.logger.info(f"  Last modified: {file_date.iso}")

                backup_file: str = self.LIST_FILE.replace(
                    '.txt', file_date.isot[:16].replace('-', '').replace(':', ''))
                os.system(f'cp {self.LIST_FILE} {backup_file}')

    def new_label_urls(self) -> None:
        """Iterator of new labels."""
        line_count: int = 0
        calibrated_count: int = 0
        processed_count: int = 0
        with open(self.LIST_FILE, "r") as inf:
            line: str
            for line in inf:
                line_count += 1
                if re.match(".*data_calibrated/.*\.xml\n$", line):
                    if "collection" in line:
                        continue
                    calibrated_count += 1
                    path: str = line.strip()
                    path = path[line.find("gbo.ast.catalina.survey"):]
                    processed: Union[None, SyncStatus] = (
                        self.db.session.query(SyncStatus)
                        .filter(SyncStatus.path == path)
                        .first()
                    )

                    if processed is None:
                        processed_count += 1
                        yield self.ARCHIVE_PREFIX + path

        self.logger.info("Processed:\n"
                         f"    {line_count} lines\n"
                         f"    {calibrated_count} calibrated data labels\n"
                         f"    {processed_count} new files\n")

    def process_label(self, url: str) -> Observation:
        """Download label, parse, and create database object."""

        label: PDS4Label = pds4_read(url, lazy_load=True, quiet=True).label
        lid: str = label.find("Identification_Area/logical_identifier").text

        tel: str = lid.split(":")[5][:3].upper()
        if tel in CatalinaBigelow._telescopes:
            obs = CatalinaBigelow()
        elif tel in CatalinaLemmon._telescopes:
            obs = CatalinaLemmon()
        elif tel in CatalinaKittPeak._telescopes:
            obs = CatalinaKittPeak()
        else:
            raise ValueError(f"Unknown telescope {tel}")

        obs.product_id = lid
        obs.mjd_start = Time(
            label.find("Observation_Area/Time_Coordinates/start_date_time").text
        ).mjd
        obs.mjd_stop = Time(
            label.find("Observation_Area/Time_Coordinates/stop_date_time").text
        ).mjd
        obs.exposure = round((obs.mjd_stop - obs.mjd_start) * 86400, 3)

        survey: ElementTree = label.find(".//survey:Survey")
        ra: List[float] = []
        dec: List[float] = []
        corner: str
        for corner in ("Top Left", "Top Right", "Bottom Right", "Bottom Left"):
            coordinate: ElementTree = survey.find(
                "survey:Image_Corners"
                f"/survey:Corner_Position[survey:corner_identification='{corner}']"
                "/survey:Coordinate"
            )
            ra.append(float(coordinate.find("survey:right_ascension").text))
            dec.append(float(coordinate.find("survey:declination").text))
        obs.set_fov(ra, dec)

        maglimit: Union[ElementTree, None] = survey.find(
            "survey:Limiting_Magnitudes"
            "/survey:Percentage_Limit[survey:Percentage_Limit='50']"
            "/survey:limiting_magnitude"
        )
        if maglimit is not None:
            obs.maglimit = float(maglimit.text)

        return obs

    def sync(self) -> None:
        """Sync local database with files available at PDS."""

        self.download_list_file()

        observations: List[Observation] = []
        tri: ProgressTriangle = ProgressTriangle(1, logger=self.logger, base=2)
        url: str
        msg: str
        failed: int = 0
        for url in self.new_label_urls():
            try:
                observations.append(self.process_label(url))
                msg = "added"
            except ValueError as e:
                failed += 1
                msg = str(e)
            except:
                self.logger.error(
                    "A fatal error occurred processing %s", url, exc_info=True
                )
                raise

            self.logger.debug("%s: %s", url, msg)
            tri.update()

            status = SyncStatus(
                path=url[len(self.ARCHIVE_PREFIX):],
                date=Time.now().iso,
                status=msg
            )
            self.db.session.add(status)

            if len(observations) >= 10000:
                self.add_observations(observations)
                observations = []

        tri.log()

        if len(observations) > 0:
            self.add_observations(observations)

        if failed > 0:
            self.logger.warning("Failed processing %d files", failed)
