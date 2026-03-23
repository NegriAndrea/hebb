#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import logging
logger = logging.getLogger(__name__)

def loggerN(msg):
    logger.info(msg)

def loggerH(msg):
    logger.info('')
    logger.info(msg)
