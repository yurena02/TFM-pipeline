#!/usr/bin/env python
# coding: utf-8


import sys
import yaml


def get_nested_value(data, keys):
    """
    Accede a claves anidadas:
    filtering.vfdb.min_identity
    """

    for key in keys:
        data = data[key]

    return data


def main():

    if len(sys.argv) != 3:
        print(
            "Uso: load_config.py <config.yaml> <clave>"
        )
        sys.exit(1)

    config_file = sys.argv[1]
    query = sys.argv[2]

    with open(config_file) as f:
        config = yaml.safe_load(f)

    keys = query.split(".")

    value = get_nested_value(config, keys)

    print(value)

if __name__ == "__main__":
    main()
