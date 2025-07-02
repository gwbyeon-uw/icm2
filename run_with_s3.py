import argparse
import boto3
import json
import logging
import os
from pathlib import Path
import parser
import subprocess
import tempfile
from urllib.parse import urlparse

logging.basicConfig(format="%(asctime)-15s: %(message)s", level=logging.INFO)


def main():
    parser = argparse.ArgumentParser(description="Runs the M2Seq pipeline")
    parser.add_argument(
        "--target_descriptor", type=str, default=os.environ.get("TARGET_DESCRIPTOR")
    )
    args = parser.parse_args()

    if not args.target_descriptor:
        raise Exception("Missing target descriptor")

    s3 = boto3.client("s3")

    with tempfile.TemporaryDirectory() as tmp:
        with open(get_file_path(args.target_descriptor, tmp, s3), "r") as f:
            descriptor = json.loads(f.read())

        out_dir = Path(tmp, "out")
        os.mkdir(out_dir)
        output_prefix = "result"
        index = Path(out_dir, "index.fa")

        with open(index, "w") as f:
            for sample in descriptor["samples"]:
                name = sample["name"].replace(" ", "_")
                f.write(f">{name}\n")
                f.write(f'{sample["sequence"]}\n')

        read_1 = get_file_path(descriptor["read_1"], tmp, s3)
        read_2 = get_file_path(descriptor["read_2"], tmp, s3)

        logging.info(f"Running m2seq")
        script = Path(os.path.dirname(os.path.realpath(__file__)), "p0.sh")
        result = subprocess.run(
            [
                script,
                read_1,
                read_2,
                out_dir,
                output_prefix,
                index,
            ],
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise Exception(f"Running m2seq failed ({result.returncode}): {result.stderr}")

        if args.target_descriptor.startswith("s3://"):
            logging.info(f"Uploading results to S3")
            destination = os.path.dirname(args.target_descriptor)
            result = subprocess.run(
                [
                    "aws",
                    "s3",
                    "cp",
                    out_dir,
                    destination,
                    "--exclude",
                    "*",
                    "--include",
                    f"{output_prefix}*",
                    "--recursive",
                ],
                capture_output=True,
                text=True,
            )

            if result.returncode != 0:
                raise Exception(f"Copying to S3 failed ({result.returncode}): {result.stderr}")


def get_file_path(uri, tmp_dir: Path, s3) -> Path:
    parsed = urlparse(uri)
    if parsed.scheme == "file":
        return Path(parsed.path)
    elif parsed.scheme == "s3":
        logging.info(f"Downloading {uri} from S3")
        filename = Path(tmp_dir, os.path.basename(parsed.path))
        s3.download_file(parsed.netloc, parsed.path[1:], str(filename))
        return filename
    else:
        raise Exception(f"File {uri} is not on a known filesystem")


if __name__ == "__main__":
    main()
