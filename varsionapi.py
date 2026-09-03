from dataclasses import dataclass
from typing import Any

import requests


API_URL = "http://192.168.102.149/daq_settings"


@dataclass
class EthernetConfig:
    dhcp: bool
    ip_cfg: str
    gateway: str
    netmask: str
    dns: str


@dataclass
class DaqStatus:
    device_pid: str
    board_position: int
    daq_ip: str
    fpga_temp: float
    daq_pid: str
    daq_version: str
    firmware_default: bool
    firmware_name: str
    firmware_version: str
    uid: str
    ethernet: EthernetConfig

    @classmethod
    def from_api_response(cls, response: dict[str, Any]) -> "DaqStatus":
        if not response.get("Result", False):
            raise RuntimeError(f"DAQ returned an error: {response}")

        if response.get("Response") != "get_daq_status":
            raise ValueError(
                f"Unexpected response type: {response.get('Response')}"
            )

        payload = response["payload"]
        ethernet = payload["ethernet"]

        return cls(
            device_pid=payload["device_pid"],
            board_position=int(payload["board_position"]),
            daq_ip=payload["daq_ip"],
            fpga_temp=float(payload["fpga_temp"]),
            daq_pid=payload["daq_pid"],
            daq_version=payload["daq_ver"],
            firmware_default=payload["daq_fw_default"] == "1",
            firmware_name=payload["daq_fw_name"],
            firmware_version=payload["daq_fw_ver"],
            uid=payload["uid"],
            ethernet=EthernetConfig(
                dhcp=bool(ethernet["dhcp"]),
                ip_cfg=ethernet["ip_cfg"],
                gateway=ethernet["gw"],
                netmask=ethernet["nm"],
                dns=ethernet["dns"],
            ),
        )


def get_daq_status() -> DaqStatus:
    request_payload = {
        "command": "get_daq_status"
    }

    response = requests.post(
        API_URL,
        json=request_payload,
        timeout=5,
    )
    response.raise_for_status()

    return DaqStatus.from_api_response(response.json())


def main() -> None:
    try:
        status = get_daq_status()

        print(f"Firmware name:    {status.firmware_name}")
        print(f"Firmware version: {status.firmware_version}")

    except requests.RequestException as exc:
        print(f"Communication error: {exc}")
    except (KeyError, TypeError, ValueError, RuntimeError) as exc:
        print(f"Invalid DAQ response: {exc}")


if __name__ == "__main__":
    main()