#!/usr/bin/env python
# -*- coding: utf-8 -*-
import glob
import os
from types import MethodType

import pandas as pd

import pycif.utils.check.errclass as error
from pycif.utils.check import info
from .baseclass import Plugin


class ObsParser(Plugin):
    """Class for handling time series parsing from different data providers
    and data file formats.

    """

    def initiate_template(self):
        module = super(ObsParser, self).initiate(plg_type="obsparser")

        # Replacing auxiliary functions
        if hasattr(module, "do_parse"):
            self.do_parse = MethodType(module.do_parse, self)

    @classmethod
    def get_parser(cls, plg):
        """Get the correct Parser for a provider and file_format_id

        Args:
            provider (str):  provider of the input file
            file_format_id (str): name of the type of file with a given format

        Returns:
            Parser: Parser for provider and file_format_id
        """

        return cls.load_registered(
            plg.provider, plg.format, "obsparser", plg_orig=plg
        )

        # plgtmp = cls.get_registered(plg.provider, plg.format, 'obsparser')
        #
        # for attr in plg.attributes:
        #     setattr(plgtmp, attr, getattr(plg, attr))
        #
        # # Adding plugin attribute
        # plgtmp.plugin = cls.from_dict(
        #     {'name': plg.provider, 'version':plg.format, 'type': 'obsparser'}
        # )
        #
        # # Creating an ObsParser instance and initializing it
        # plgtmp = cls(plg_orig=plgtmp)
        # plgtmp.initiate_template()
        #
        # return plgtmp

    @classmethod
    def register_parser(cls, provider, file_format_id, parse_module, **kwargs):
        """Register a parsing function for provider and format with default
        options

        Args:
            provider (str):  provider of the input file
            file_format_id (str): name of the type of file with a given format
            parse_module (Module):
                    returns file content
                    as pandas.DataFrame df[obssite_id, parameter]
            **kwargs: default options for parse_function

        Notes:
            The parse_function signature is the same as
            the :py:func:`Parser.parse_file`
        """
        super(ObsParser, cls).register_plugin(
            provider, file_format_id, parse_module, plugin_type="obsparser"
        )

    def parse_file(self, obs_file, **kwargs):
        """This function does the parsing (and post processing if necessary).

        Args:
            obs_file (str): path to input file

        Keyword Args:
            encoding (str): Encoding of input files
            freq (str): frequency after resampling
                        see `Offset Aliases`_ for valid strings
            src_freq (str):
                        explicit setting of the frequency in the input file
                        shouldn't be necessary

        Returns:
            pandas.DataFrame: renamed, shifted, resampled
            Dataframe df[obssite_id, parameter] with t as index
        """

        df = self.do_parse(obs_file, **kwargs)

        # Removing rows with only NaNs
        df = df.dropna(axis=1, how="all")

        # Checking that the returned dataframe as all required columns
        if self.check_df(df, **kwargs):
            return df

    def parse_multiple_files(self, **kwargs):
        """Parses multiple files specified by a glob pattern and stores the
        content into a datastore

        Args:
            provider_name (str):  provider of the input file
            file_format_id (str): name of the type of file with a given format
            glob_pattern (str): glob pattern: /**/ for recursive  matching
                                in subdirectories

        Keyword Args:
            encoding (str): Encoding of input files
            freq (str): frequency after re-sampling
                        see `Offset Aliases`_ for valid strings
            src_freq (str): explicit setting of the frequency in the input file
                            shouldn't be necessary

        Notes:
            - Additional kwargs for a parser are possible or even required.
              See the respective documentation

        Returns:
            dict: {obs_file} = df[obssite_id, parameter]
        """

        # parser = cls.get_parser(provider_name, file_format_id)

        dfs = {}

        info("Reading files in " + self.dir_obs)

        for obs_file in sorted(glob.glob(self.dir_obs + "*")):
            try:
                dfs[os.path.basename(obs_file)] = self.parse_file(
                    obs_file, **kwargs
                )

            except error.PluginError as e:
                info(
                    "{} was not loaded for the following reason".format(
                        obs_file
                    )
                )
                info(e.message)

        if dfs != {}:
            return pd.concat(list(dfs.values()))
        else:
            return pd.DataFrame({})

    @staticmethod
    def check_df(df, **kwargs):

        reqcols = ["station", "network", "parameter", "duration", "obserror"]
        cols = df.columns

        for c in reqcols:
            if c not in cols:
                towrite = """
                    {} was not returned in the dataframe
                    Please check your parser definition
                    """.format(
                    c
                )

                raise error.PluginError(towrite)

        return True
