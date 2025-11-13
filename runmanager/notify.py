#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__date__ = '27 Aug. 2025'

import requests
import logging

logger = logging.getLogger(__name__)

def send_to_discord(webhook_url: str, message: str):
    """
    Sends a message to the specified Discord webhook URL.
    
    Args:
        webhook_url (str): The Discord webhook URL.
        message (str): The message content to send.
    """
    if not webhook_url:
        # Silently ignore if webhook_url is not provided (e.g., in viewer mode).
        return
    try:
        payload = {"content": message}
        response = requests.post(webhook_url, json=payload, timeout=10)
        response.raise_for_status()  # Raise an exception for bad status codes
        logger.info('Successfully sent a notification to Discord.')
    except requests.exceptions.RequestException as e:
        logger.error(f'Failed to send notification to Discord: {e}')
