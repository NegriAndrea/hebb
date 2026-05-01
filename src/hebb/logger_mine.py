#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[
        # logging.FileHandler("hebb_run.log"),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger(__name__)

def loggerN(msg):
    logger.info(msg)

def loggerH(msg):
    logger.info('')
    logger.info(msg)
