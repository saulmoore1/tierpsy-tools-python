#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 15:52:14 2020

@author: em812
"""

import pandas as pd
from pathlib import Path
from warnings import warn
import pdb

def read_hydra_metadata(
        feat_file, fname_file, meta_file,
        feat_id_cols = ['file_id', 'n_skeletons', 'well_name', 'is_good_well'],
        add_bluelight=True,
        bluelight_labels=['prestim', 'bluelight', 'poststim']):

    """
    Creates matching features and metadata dfs from the .csv files of a hydra
    screening (assuming the standardized format of tierpsy and hydra metadata
    from tierpsytools).

    Parameters
    ----------
    feat_file : file path to tierpsy features summaries file
        File must have a file_id and well_id column.
    fname : file path to tierpsy filenames summaries file.
        File must have a file_id and filename column.
    meta : file path to metadata file
    feat_id_cols : list of strings, optional
        The columns in the feat dataframe that are not features.
        The default is ['file_id', 'well_name', 'is_good_well'].
    add_bluelight : bool, optional
        Add a metadata column that specifies the bluelight condition for each row.
        The default is True.
    bluelight_labels : list, optional
        The names of the bluelight conditions as they appear in the file names.
        Only used if add_bluelight is True.

    Returns
    -------
    feat: dataframe shape = (n_wells_with_features*n_bluelight_conditions, n_features)
        Features dataframe containing only feature columns..
    newmeta : dataframe shape = (n_wells_with_features*n_bluelight_conditions, n_meta_cols)
        Metadata dataframe matching the returned feat dataframe row-by-row.

    """
    feat, fname, meta = _read_files(feat_file, fname_file, meta_file, comment='#')

    feat,meta = build_matching_feat_meta(
        feat, fname, meta, feat_id_cols = feat_id_cols,
        add_bluelight=add_bluelight, bluelight_labels=bluelight_labels)

    return feat, meta


def build_matching_feat_meta(
        feat, fname, meta,
        feat_id_cols = ['file_id', 'n_skeletons', 'well_name', 'is_good_well'],
        add_bluelight=True,
        bluelight_labels=['prestim', 'bluelight', 'poststim'],
        merge_how='left'):
    """
    Creates matching features and metadata dfs from the .csv files of a hydra
    screening (assuming the standardized format of tierpsy and hydra metadata
    from tierpsytools).

    Parameters
    ----------
    feat : dataframe size =
           (n_wells_with_features*n_bluelight_conditions,
           n_features + len(feat_id_cols))
        The tierpsy features_summaries read into a dataframe.
        Must have a file_id and well_id column.
    fname : dataframe size = (n_files, n_cols)
        The tierpsy filenames_summaries read into a dataframe. Must have a
        file_id and filename column.
    meta : dataframe size = (n_wells_screened*n_bluelight_conditions, n_meta_cols)
        The experiment full metadata read into a dataframe.
    feat_id_cols : list of strings, optional
        The columns in the feat dataframe that are not features.
        The default is ['file_id', 'well_name', 'is_good_well'].
    add_bluelight : bool, optional
        Add a metadata column that specifies the bluelight condition for each row.
        The default is True.
    bluelight_labels : list, optional
        The names of the bluelight conditions as they appear in the file names.
        Only used if add_bluelight is True.
    merge_how: string, optional
        Not recommended to change it.

    Returns
    -------
    feat: dataframe shape = (n_wells_with_features*n_bluelight_conditions, n_features)
        Features dataframe containing only feature columns..
    newmeta : dataframe shape = (n_wells_with_features*n_bluelight_conditions, n_meta_cols)
        Metadata dataframe matching the returned feat dataframe row-by-row.

    """
    if 'filename' in fname:
        filename = 'filename'
    elif 'file_name' in fname:
        filename = 'file_name'
    else:
        raise ValueError('The filenames dataframe needs to have a filename column.')

    fname['imgstore_name'] = fname[filename].apply(
        lambda x: _imgstore_name_from_filename(x,path_levels=[-3,-1]))

    newmeta = feat[feat_id_cols]

    newmeta.insert(
        0, 'imgstore_name', newmeta['file_id'].map(
            dict(fname[['file_id', 'imgstore_name']].values))
        )
    newmeta.insert(
        0, 'featuresN_filename', newmeta['file_id'].map(
            dict(fname[['file_id', filename]].values))
        )

    newmeta = pd.merge(
        newmeta, meta, on=['imgstore_name','well_name'], how=merge_how
        )

    if add_bluelight:
        newmeta = add_bluelight_label(
            newmeta, labels=bluelight_labels)

    if merge_how=='left':
        assert newmeta.shape[0] == feat.shape[0]
    elif merge_how=='right':
        if not newmeta.shape[0] == meta.shape[0]:
            breakpoint()

    _cols = feat_id_cols+['imgstore_name', 'featuresN_filename']
    if add_bluelight:
        _cols += ['bluelight']
    all_nans = newmeta[newmeta.columns.difference(_cols)].isna().all(axis=1)
    if  all_nans.sum() > 1:
        warn(
            'There are wells in the feature summaries that do not have '+
            'matching metadata information. These wells will be discarded.')
        newmeta = newmeta[~all_nans]
        feat = feat[~all_nans]
        newmeta = newmeta.reset_index(drop=True)
        feat = feat.reset_index(drop=True)

    return feat[feat.columns.difference(feat_id_cols)], newmeta

def align_bluelight_conditions(
        feat, meta,
        how = 'outer',
        return_separate_feat_dfs = False,
        bluelight_specific_meta_cols = ['imgstore_name', 'file_id', 'bluelight', 'n_skeletons'],
        merge_on_cols = ['date_yyyymmdd','imaging_plate_id','well_name']
        ):
    """
    Concatenates the features from the three bluelight conditions for each well
    and creates matching features and metadata dataframes for the concatenated
    feature matrix.

    !: is bad well not considered bluelight specific!

    Parameters
    ----------
    feat: dataframe shape = (n_wells_with_features*n_bluelight_conditions, n_features) or None
        Features dataframe containing only feature columns.
        It can be set to None if we only want to align the metadata dataframe.
    meta : dataframe shape = (n_wells_with_features*n_bluelight_conditions, n_meta_cols)
        Metadata dataframe matching the feat dataframe row-by-row.
    how : string, optional
        Specifies merge type. The default is 'outer'.
    return_separate_feat_dfs : bool, optional
        If True, three different dataframes will be returned, one for each
        bluelight condition. The default is False.
    bluelight_specific_meta_cols : list of strings, optional
        The metadata columns that have different values for each bluelight condition.
        The default is ['imgstore_name', 'file_id', 'bluelight'].
    merge_on_cols : list of strings, optional
        The metadata columns that together point to a unique well.
        The default is ['date_yyyymmdd','imaging_plate_id','well_name'].

    Returns
    -------
    IF return_separate_feat_dfs:
        feat: dataframe shape = (n_wells, n_features * n_bluelight_conditions)
            The features at every bluelight condition for each well.
        meta: dataframe shape = (n_wells, n_meta_cols)
            The metadata tht match the returned feat dataframe row-by-row

    ELSE:
        (feat_prestim, feat_bluelight, feat_poststim): tuple of dataframes
                each one with shape = (n_wells, n_features)
            The features at every bluelight condition for each well.
        meta: dataframe shape = (n_wells, n_meta_cols)
            The metadata that match the returned feat dataframes row-by-row.

    """

    bluelight_conditions = ['prestim', 'bluelight', 'poststim']

    # The metadata columns that are the same for all bluelight conditions
    shared_meta_cols = meta.columns.difference(
        bluelight_specific_meta_cols+merge_on_cols).to_list()

    if feat is not None:
        # Store the feature names
        feat_names = feat.columns.to_list()

        # Concatenate the metadata with the features for a uniform merge based on
        # metadata
        meta = pd.concat([meta, feat], axis=1)
    else:
        feat_names = []

    # Merge based on the merge_on_cols columns
    meta = meta.set_index(merge_on_cols)

    meta = meta[meta['bluelight']=='prestim'].join(
            meta[meta['bluelight']=='bluelight'],
            how=how,
            lsuffix='',
            rsuffix='_bluelight').join(
                meta[meta['bluelight']=='poststim'],
                how=how,
                lsuffix='',
                rsuffix='_poststim')

    # Add the prestim suffix to the bluelight-specific columns
    meta = meta.rename(
            columns={ft:ft+'_prestim' for ft in bluelight_specific_meta_cols+feat_names})

    if feat is not None:
        # Derive the feature column names in the merged dataframe and use them to
        # split again the features from the metadata
        feat_cols = ['_'.join([ft,blue])
                     for blue in bluelight_conditions
                     for ft in feat_names]

        meta.reset_index(drop=False, inplace=True)
        feat = meta[feat_cols]
        meta = meta[meta.columns.difference(feat_cols)]


    # The shared meta columns might have nan values in some of the conditions
    # because of the outer merge. I want to keep only one of these columns for
    # all the conditions. So I will replace the nans in the prestim version
    # with the non-nan values from the other conditions and then drop the
    # columns of the other conditions.
    for col in shared_meta_cols:
        for blue in bluelight_conditions[1:]:
            meta.loc[meta[col].isna(), col] = \
                meta.loc[meta[col].isna(), '_'.join([col, blue])]
            meta = meta.drop(columns=['_'.join([col, blue])])

    # Return
    if return_separate_feat_dfs:
        prestim = [ft for ft in feat.columns if '_prestim' in ft]
        bluelight = [ft for ft in feat.columns if '_bluelight' in ft]
        poststim = [ft for ft in feat.columns if '_poststim' in ft]
        return (feat[prestim], feat[bluelight], feat[poststim]), meta
    else:
        return feat, meta


# %% helper functions
def add_bluelight_label(
        meta, filename_column='imgstore_name', split_string='_',
        labels=['prestim', 'bluelight', 'poststim']
        ):

    from numpy import nan

    if 'bluelight' in meta:
        return meta

    def _bluelight_from_filename(fname):
        split = Path(fname).stem.split(split_string)
        isin = [x in split for x in labels]
        label = [l for (l, i) in zip(labels, isin) if i]
        if len(label) == 1:
            return label[0]
        else:
            return nan

    bluelight = meta[filename_column].apply(_bluelight_from_filename)
    meta.insert(0, 'bluelight', bluelight)
    return meta

def _imgstore_name_from_filename(filename, path_levels=[-3,-1]):

    imgstore_name = '/'.join(Path(filename).parts[path_levels[0]:path_levels[1]])
    return imgstore_name

def _read_files(feat_file, fname_file, metadata_file, comment='#'):
    feat = pd.read_csv(feat_file, comment='#')
    fname = pd.read_csv(fname_file, comment='#')
    meta = pd.read_csv(metadata_file, index_col=None)

    return feat, fname, meta