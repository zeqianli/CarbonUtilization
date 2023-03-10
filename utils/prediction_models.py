""" Utilities for the carbon utilization project. 

@author: Zeqian Li
@contact: zeqianli.chicago@gmail.com
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import itertools
import random
import traceback
import scipy.cluster.hierarchy as hierarchy
from functools import reduce
from tqdm import tqdm
from multiprocessing import Pool
from ete3 import Tree
from Bio import Phylo
from io import StringIO
from scipy.stats import ttest_ind
from sklearn.base import BaseEstimator
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, accuracy_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from statannotations.Annotator import Annotator
import sklearn.metrics as metrics 
from scipy.stats import entropy
from statsmodels.stats.multitest import multipletests


plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize'] = (7, 5)
plt.rcParams['axes.axisbelow'] = True

pd.set_option('display.max_columns', None)

# Carbon sources in Zeqian dataset.
CARBONS = ['Arabinose', 'Butyrate','Deoxyribose','Glucuronic acid', 'Glycerol', 'Mannitol', 'Mannose', 'Melibiose', 'Propionate', 'Raffinose']

# default colors
COLORS=sns.color_palette("colorblind").as_hex()
COLORS.pop(7) # remove grey

# ============================================================
# Prediction models
# ============================================================

class BinaryGrowthClassifier(BaseEstimator):
    """ Base class for all binary growth classifiers.

    Child classes must implement four methods: __init__, fit, predict, and get_params.
    """

    def __init__(self, **kwargs):
        return

    def fit(self, X, y):
        """ Fit the model. Must be implemented by child classes.

        Parameters
        ----------
        X: an array of shape (n_samples, n_features).
        y: a binary array of shape (n_samples,). Must not contrain any NaN values.
        """
        raise NotImplementedError

    def predict(self, X):
        """ Made predictions. Must be implemented by child classes.

        Returns
        -------
        y: a binary array of shape (n_samples,). 
        """
        raise NotImplementedError

    def score(self, X, y):
        """ Compute the prediction score. Children class don't need to implement this method.

        Parameters
        ----------
        X: an array of shape (n_samples, n_features).
        y: a binary array of shape (n_samples,). Must not contrain any NaN values.

        Returns
        -------
        score: a dictionary of these scores:
            - true_negative, false_positive, false_negative, true_positive
            - accuracy: ratio of correct prediction. Never NaN.
            - precision: may be NaN. 
            - sensitivity: may be NaN.
            - specificity: may be NaN.
            - f1_score: may be NaN.
            - balanced_accuracy: may be NaN.
        """
        prediction = self.predict(X)
        cm = confusion_matrix(y, prediction, labels=[0, 1])
        score = {'true_negative': cm[0, 0], 'false_positive': cm[0, 1],
                 'false_negative': cm[1, 0], 'true_positive': cm[1, 1]}

        score['accuracy'] = accuracy_score(
            y, prediction)  # TODO: deal with zeros

        if (score['true_positive']+score['false_positive']) > 0:
            score['precision'] = score['true_positive'] / \
                (score['true_positive']+score['false_positive'])
        else:
            score['precision'] = np.nan

        if (score['true_positive']+score['false_negative']) > 0:
            score['sensitivity'] = score['true_positive'] / \
                (score['true_positive']+score['false_negative'])
        else:
            score['sensitivity'] = np.nan

        if (score['true_negative']+score['false_positive']) > 0:
            score['specificity'] = score['true_negative'] / \
                (score['true_negative']+score['false_positive'])
        else:
            score['specificity'] = np.nan

        if (not np.isnan(score['precision'])) and (not np.isnan(score['sensitivity'])) and (score['precision']+score['sensitivity']) > 0:
            score['f1_score'] = 2*score['precision'] * \
                score['sensitivity']/(score['precision']+score['sensitivity'])
        else:
            score['f1_score'] = np.nan
        if (not np.isnan(score['sensitivity'])) and (not np.isnan(score['specificity'])):
            score['balanced_accuracy'] = (
                score['sensitivity']+score['specificity'])/2
        else:
            score['balanced_accuracy'] = np.nan
        # score['recall']=recall_score(y,prediction)
        # score['precision']=precision_score(y,prediction)
        # score['f1_score']=f1_score(y,prediction)
        return score

    def get_params(self, deep=True):
        """ Return a dictionary of parameters to be recorded. Must be implemented by child classes."""
        raise NotImplementedError


class RF(BinaryGrowthClassifier):
    """ Random forest classifier.

    Parameters: see sklearn.ensemble.RandomForestClassifier. These parameters are common:
        - n_estimators: number of trees in the forest. Default: 100.
        - max_depth: maximum depth of the tree. Default: None.
        - max_features: number of features to consider when looking for the best split. Default: 'sqrt'.
    """

    def __init__(self, **kwargs):
        self.model = RandomForestClassifier(**kwargs)
        return

    def fit(self, X, y):
        """ Fit the random forest

        Parameters
        ----------
        X: an array of shape (n_samples, n_features).
        y: a binary array of shape (n_samples,). Must not contrain any NaN values.
        """
        self.model.fit(X, y)
        return self

    def predict(self, X):
        """ Made predictions. Must be implemented by child classes.

        Returns
        -------
        y: a binary array of shape (n_samples,). 
        """
        return self.model.predict(X)

    def get_params(self, deep=True):
        """ Return these parameters:

        - feature_importances: an array of feature importances of shape (n_features,).
        """
        return {'feature_importances': self.model.feature_importances_}


class DecisionTree(BinaryGrowthClassifier):
    def __init__(self, **kwargs):
        self.model = DecisionTreeClassifier(**kwargs)

    def fit(self, X, y):
        self.model.fit(X, y)
        return self

    def predict(self, X):
        return self.model.predict(X)

    def get_params(self, deep=True):
        return {'feature_importances': self.model.feature_importances_,
                'max_features': self.model.max_features_,
                'max_depth': self.model.max_depth}


class FBAClassifier(BinaryGrowthClassifier):

    bigg_id_map = pd.DataFrame([
        ['Arabinose', 'EX_arab__L_e'],
        ['Deoxyribose', 'EX_drib_e'],
        ['Glucuronic acid', 'EX_glcur_e'],
        ['Glycerol', 'EX_glyc_e'],
        ['Mannitol', 'EX_mnl_e'],
        ['Mannose', 'EX_man_e'],
        ['Melibiose', 'EX_melib_e'],
        ['Raffinose', 'EX_raffin_e'],
        ['Butyrate', 'EX_but_e'],
        ['Propionate', 'EX_ppa_e']
    ], columns=['carbon_name', 'bigg_id']).set_index('carbon_name')

    def __init__(self, growth, c=None, growth_threshold=0.01):
        if isinstance(growth, str):
            # .drop(index=['EX_mobd_e']) # I don't know why there is an extra row at the end
            growth = pd.read_csv(growth, index_col=0)
            growth = growth.loc[FBAClassifier.bigg_id_map.at[c, 'bigg_id'], :]
        self.growth = growth
        self.growth_threshold = growth_threshold

    def fit(self, X, y):
        return self

    def predict(self, X):
        return (self.growth > self.growth_threshold).loc[X.index.values].astype(int)

    def get_params(self, deep=True):
        return {'growth_threshold': self.growth_threshold}


class NullClassifier(BinaryGrowthClassifier):
    """ Base class for null models. 

    Child classes must implement four methods: __init__, fit, predict, get_params.
    """

    def get_params(self, deep=True):
        """ Return a dictionary of model parameters. """
        return {'is_null_model': True}


class BernoulliNull(NullClassifier):
    """ Bernoulli null model. Generates random predictions with the same mean as the training data. """

    def __init__(self):
        self.p = None

    def fit(self, X, y):
        self.p = y.mean()
        return self

    def predict(self, X):
        return np.random.binomial(1, self.p, size=X.shape[0])

    def get_params(self, deep=True):
        out = super().get_params(deep=deep)
        out.update({'null_model': 'bernoulli'})
        return out

class IdentityNull(NullClassifier):
    """ Identity null model. Predicts the same value for all samples. """

    def __init__(self):
        self.pred = None

    def fit(self, X, y):
        self.pred = int(y.mean()>=0.5)
        return self

    def predict(self, X):
        return self.pred * np.ones(X.shape[0])

    def get_params(self, deep=True):
        out = super().get_params(deep=deep)
        out.update({'null_model': 'bernoulli', 'identity_null_pred': self.pred})
        return out



class RandomGrowthNull(NullClassifier):
    """ Fit models on random growth data. 
    
    Parameters
    ----------
    Model: a BinaryGrowthClassifier class.

    replace: bool, default=True
        When generating random growth data, whether to sample with replacement.
    
    **kwargs: arguments passed to Model.
    """
    def __init__(self, Model, replace=True, **kwargs):
        self.model = Model(**kwargs)
        self.replace = replace

    def fit(self, X, y):
        y_shuffled = pd.Series(np.random.choice(
            y.values, size=len(y), replace=self.replace), index=y.index)
        self.model.fit(X, y_shuffled)
        return self

    def predict(self, X):
        return self.model.predict(X)

    def get_params(self, deep=True):
        out = super().get_params(deep=deep)
        out.update({'null_model': 'random_growth', 'replace': self.replace})
        out.update(self.model.get_params())
        return out


class PretrainedModel(BinaryGrowthClassifier):
    def __init__(self, model):
        self.model = model
        return

    def fit(self, X, y):
        return self

    def predict(self, X):
        return self.model.predict(X)

    def get_params(self, deep=True):
        out = {'pre-trained': True}
        out.update(self.model.get_params())
        return out


class NearestNeighbor(BinaryGrowthClassifier):
    """ K-nearest neighbor classifier. 
    
    Parameters
    ----------
    n_neighbors: int, default=1
        Number of neighbors to use by default for kneighbors queries.
    
    **kwargs: arguments passed to KNeighborsClassifier.
    """
    def __init__(self, n_neighbors=1, **kwargs):
        self.model = KNeighborsClassifier(n_neighbors=n_neighbors, **kwargs)

    def fit(self, X, y):
        self.model.fit(X, y)
        return self

    def predict(self, X):
        return self.model.predict(X)

    def get_params(self, deep=True):
        return {'n_neighbors': self.model.n_neighbors,
                'weights': self.model.weights,
                'metric': self.model.effective_metric_}


class LassoLogisticRegression(BinaryGrowthClassifier):
    """ Lasso logistic regression classifier.

    Parameters
    ----------
    penalty: str, default='l1'
        Used to specify the norm used in the penalization. The 'newton-cg', 'sag' and 'lbfgs' solvers support only l2 penalties.

    C: float, default=1.0
        Inverse of regularization strength; must be a positive float. Like in support vector machines, smaller values specify stronger regularization.

    solver: str, default='liblinear'
        Algorithm to use in the optimization problem.
    
    **kwargs: arguments passed to LogisticRegression.
    """ 
    def __init__(self, penalty='l1', C=1.0, solver='liblinear', **kwargs):
        self.model = LogisticRegression(
            penalty=penalty, C=C, solver=solver, **kwargs)

    def fit(self, X, y):
        self.model.fit(X, y)
        return self

    def predict(self, X):
        return self.model.predict(X)

    def get_params(self, deep=True):
        return {'coefficients': self.model.coef_, 'intercept': self.model.intercept_, 'C': self.model.C}


class GreedyFeatureSelection(BinaryGrowthClassifier):
    """ Greedy feature selection classifier. 

    Parameters
    ----------
    Model: a BinaryGrowthClassifier class.

    model_params: dict, default=None
        Parameters to pass to Model.

    tree: ete3.Tree, default=None
        Required if split_method=='ooc'. 
    
    n_max_features: int, default=5
        Maximum number of selected features.

    n_feature_subsample: int, default=None
        If not None, subsample features in each searching step to accelerate the selection process. 
    
    p: Pool, default=None
        Required if threads>1. 
    
    save_meta_model: bool, default=False
        If True, save the meta-models objects in the final data. Be careful about memory usage. 
    
    multithreading_batch: int or None, default=None 



    """

    def __init__(self,
                 Model, 
                 model_params=None, 
                 tree=None,
                 n_max_features=5, 
                 n_feature_subsample=None, 
                 # improvement_cutoff=0.05,
                 n_meta_split=10, 
                 split_method='ooc', 
                 splitter_params=None,
                 pick_feature_metric='accuracy',
                 threads=1, 
                 p=None, 
                 ff_results=None,
                 save_meta_models=False, 
                 multithreading_batch=None, 
                 save_meta_training_results=True,
                 verbose=True):
        self.tree = tree

        self.Model = Model
        self.model_params = model_params
        self.n_max_features = n_max_features
        self.n_feature_subsample = n_feature_subsample

        self.n_meta_split = n_meta_split
        self.split_method = split_method
        self.splitter_params = splitter_params

        self.pick_feature_metric = pick_feature_metric
        self.threads = threads
        if multithreading_batch is None:
            multithreading_batch = threads * 5
        self.multithreading_batch = multithreading_batch
        # if p is None:
        #     self._external_pool = False
        # else:
        #     self._external_pool = True
        if threads>1 and p is None:
            raise ValueError("A Pool object is required for multithreading.")
        self._p = p  # multiprocessing pool
        self.save_meta_models = save_meta_models
        # self.improvement_cutoff = improvement_cutoff

        self.verbose=verbose

        # Caching
        if ff_results is not None and os.path.exists(ff_results):
            print(f"Warning: output file {ff_results} exists. ")
            i = 0
            _spl = ff_results.split('.')
            while os.path.exists(ff_results):
                ff_results = '.'.join(_spl[:-1])+'_'+str(i)+'.'+_spl[-1]
                i += 1
            print(f"Changed to {ff_results}. ")
        self.ff_results = ff_results
        self._results = []

        self.save_meta_training_results = save_meta_training_results
        self._meta_training_data = []
        self._best_features, self._best_accuracies = None, None
        self.final_features, self.final_accuracy, self.final_model = None, None, None

    def fit(self, X, y):
        # generate splits
        splits = []
        samples = X.index.values
        if isinstance(self.split_method, type):
            splitter = self.split_method(**self.splitter_params)
        elif self.split_method == 'random':
            splitter = RandomSplitter(**self.splitter_params)
        elif self.split_method == 'ooc':
            splitter = LargeTreeTraverseOOCSplitter(
                tree=self.tree, growth_data=y, **self.splitter_params)
        # elif self.split_method == 'ooc_distmat':
        #     splitter = DistMatOOCSplitter(
        #         tree=self.tree, growth_data=y, **self.splitter_params)
        elif self.split_method == 'leave_one_out':
            splitter = LeaveOneOutSplitter(**self.splitter_params)
        for test_samples in splitter.generate_splits(samples, self.n_meta_split):
            splits.append((np.setdiff1d(samples, test_samples), test_samples))

        # run meta training

        best_features, best_accuracies = [], []
        final_features, final_accuracies = None, None

        # if self.threads > 1 and self._p is None:
        #     self._p = Pool(self.threads)
        #     print("New pool started. ")
        for n_feature in range(self.n_max_features):
            if self.n_feature_subsample is None or self.n_feature_subsample > X.shape[1]:
                iter_ = X.columns
            else:
                iter_ = np.random.choice(
                    X.columns, self.n_feature_subsample, replace=False)
            if self.verbose:
                iter_=tqdm(iter_, desc="Iterating features...")
            _data_split_batch = []
            # results=[]

            # Run meta models on all data splits and all possible features
            for i_feature, new_feature in enumerate(iter_):

                if len(best_features) > 0:
                    if new_feature in best_features[-1]:
                        continue
                    else:
                        features = best_features[-1]+[new_feature]
                else:
                    features = [new_feature]

                X_new = X[features]
                for i_split, (train_samples, test_samples) in enumerate(splits):
                    X_train = X_new.loc[train_samples]
                    y_train = y.loc[train_samples]
                    X_test = X_new.loc[test_samples]
                    y_test = y.loc[test_samples]

                    param = [self.Model(**self.model_params),
                             X_train, X_test, y_train, y_test]
                    param.append({'train_samples': train_samples, 'test_samples': test_samples,
                                 'n_feature': n_feature+1, 'features': features, 'new_feature': new_feature})
                    param.append(self.save_meta_models)

                    if self.threads <= 1:
                        self._save_new_results([_pickleable_run_model(*param)])
                    else:
                        _data_split_batch.append(param)
                        if (len(_data_split_batch) >= self.multithreading_batch) or (i_feature == len(iter_)-1 and i_split == len(splits)-1):

                            result = self._p.starmap(
                                _pickleable_run_model, _data_split_batch)
                            self._save_new_results(result)
                            _data_split_batch = []
                        else:
                            continue

            # select for best features
            results = self._load_final_results()
            self._clear_cache()
            self._meta_training_data.append(results)
            results = results.groupby('new_feature')[
                self.pick_feature_metric].mean().sort_values(ascending=False)
            best_new_feature = results.index.values[0]
            best_new_accuracy = results.values[0]

            if len(best_features) > 0:
                best_features.append(best_features[-1]+[best_new_feature])
                best_accuracies.append(best_new_accuracy)
            else:
                best_features.append([best_new_feature])
                best_accuracies.append(best_new_accuracy)

            # if len(best_accuracies) > 1 and (best_accuracies[-1]/best_accuracies[-2]) > self.improvement_cutoff:
            #     final_features, final_accuracy = best_features[-1], best_accuracies[-1]

        # if self.threads > 1 and not self._external_pool:
        #     self._p.close()
        #     print("Pool closed.")
        #     self._p = None

        if final_features is None:
            final_features, final_accuracy = best_features[-1], best_accuracies[-1]

        self._meta_training_data = pd.concat(
            self._meta_training_data, axis=0, ignore_index=True)
        self._best_features = best_features
        self._best_accuracies = best_accuracies
        self.final_features = final_features
        self.final_accuracy = final_accuracy
        self.final_model = self.train_final_model(X, y, self.final_features)

    def _save_new_results(self, results):
        if not isinstance(results, list):
            results = [results]

        if self.ff_results is None:
            self._results.extend(results)
        else:
            with open(self.ff_results, 'ab') as f:
                pickle.dump(results, f)

    def _clear_cache(self):
        if self.ff_results is None:
            self._results = []
        else:
            if os.path.exists(self.ff_results):
                os.remove(self.ff_results)

    def _load_final_results(self):
        if self.ff_results is None:
            results = self._results
        else:
            results = []
            with open(self.ff_results, 'rb') as f:
                while True:
                    try:
                        results.extend(pickle.load(f))
                    except EOFError:
                        break
        return pd.DataFrame(results)

    def train_final_model(self, X, y, final_features):
        model = self.Model(**self.model_params)
        model.fit(X[final_features], y)
        return model

    def predict(self, X):
        return self.final_model.predict(X[self.final_features])

    def get_params(self, deep=True):
        return {'final_features': self.final_features, 'final_accuracy': self.final_accuracy}


def _entropy(arr):
    if np.all(arr==0) or np.all(arr==1):
        return 0
    else:
        return entropy(arr)

def neg_conditional_entropy(arr1, arr2):
    if isinstance(arr1,list):
        arr1=np.array(arr1)
    if arr1.ndim!=1:
        arr1=(arr1 * np.array([[2**j for j in range(arr1.shape[1])] for i in range(arr1.shape[0])])).sum(axis=1) # binary coding
    vs, counts = np.unique(arr1, return_counts=True)
    ce=0

    for v, count in zip(vs, counts):
        ce+=-count/len(arr1)*_entropy(arr2[arr1==v])    
    return ce

def _pickleable_cal_nce(key, explainer, y):
    nce=neg_conditional_entropy(explainer, y)
    return key, nce

class NCEFeatureSelection(BinaryGrowthClassifier):

    def __init__(self,
                 Model, model_params, 
                 keep_top=5, max_features=5, n_feature_subsample=None, 
                 threads=1, p=None, 
                 # ff_results=None,
                 # multithreading_batch=None, 
                 save_feature_selection_data=True,
                 verbose=False):
        self.Model = Model
        self.model_params = model_params
        self.keep_top = keep_top
        self.max_features = max_features
        self.n_feature_subsample = n_feature_subsample

        self.threads = threads 
        if threads>1 and p is None:
            raise ValueError("A Pool object is required for multithreading.")
        self._p = p  # multiprocessing pool
        
        self.save_feature_selection_data  = save_feature_selection_data
        self._feature_selection_data = []

        self.verbose=verbose

        # Caching
        # if ff_results is not None and os.path.exists(ff_results):
        #     print(f"Warning: output file {ff_results} exists. ")
        #     i = 0
        #     _spl = ff_results.split('.')
        #     while os.path.exists(ff_results):
        #         ff_results = '.'.join(_spl[:-1])+'_'+str(i)+'.'+_spl[-1]
        #         i += 1
        #     print(f"Changed to {ff_results}. ")
        # self.ff_results = ff_results
        # self._results = []

        self._best_features, self._best_accuracies = None, None
        self.final_features, self.final_features_nce, self.final_model = None, None, None

    def fit(self,X,y):
        df=self.get_top_features(X, y) # feature selection data 
        if self.save_feature_selection_data:
            self._feature_selection_data = df
        self.final_features=df[df['n_features']==self.max_features].sort_values('nce')['features'].values[0]
        self.final_features_nce=np.array([df[df['features_str']==(','.join(self.final_features[:(_+1)]))]['nce'].values[0] for _ in range(self.max_features)])

        # train final model
        self.final_model=self.Model(**self.model_params)
        self.final_model.fit(X[self.final_features], y)

    def predict(self, X):
        return self.final_model.predict(X[self.final_features])

    def get_params(self, deep=True):
        return {'final_features': self.final_features, 'final_features_nce': self.final_features_nce}

    def get_top_features(self,X, y):
        if self.keep_top!=1 and isinstance(self.keep_top, float):
            self.keep_top=int(self.keep_top*X.shape[1])
        features=X.columns
        if isinstance(y,pd.Series):
            y=y.values
        past_best=None
        df=[]
        
        _iter=range(1, self.max_features+1)
        if self.verbose:
            _iter=tqdm(_iter, desc="Number of features")
        for i in _iter:
            _batch=[]
            if self.n_feature_subsample is None:
                features=X.columns
            else:
                features=np.random.choice(X.columns, self.n_feature_subsample, replace=False)
            for feature in features:
                if past_best is None:
                    _batch.append(([feature], X[feature].values, y))
                else:
                    for past_features in past_best:
                        if not feature in past_features:
                            new_features=[*past_features, feature]
                            _batch.append((new_features,X[new_features].values, y))
            if self._p is None:
                _res=[_pickleable_cal_nce(*args) for args in _batch]
            else:
                _res=self._p.starmap(_pickleable_cal_nce, _batch)
            
            _res=pd.DataFrame(_res, columns=["features", "nce"]).sort_values("nce", ascending=False)
            if isinstance(self.save_feature_selection_data, int):
                _res=_res.iloc[:self.save_feature_selection_data, :]
            past_best=_res.iloc[:self.keep_top, :]['features'].values
            df.append(_res)
        df=pd.concat(df,axis=0,ignore_index=True)
        df['n_features']=df['features'].apply(lambda x: len(x))
        df['features_str']=df['features'].apply(lambda x: ','.join(x))
        df=df.sort_values(['n_features',"nce"], ascending=[True,False])
        return df


# ============================================================
# Data set split
# ============================================================


def cal_dist_matrix_efficient(tree, df=True):
    """ Calculate pairwise distance matrix from a tree efficently using a recursive algorithm. 

    For each clade, calculate the distance matrix for each sub slade and record each of their leaf - node distance. Then add the cross-clade distances by adding the leaf - node distances of leaves from each clade.  

    Parameters:
    -----------
    tree: ete3.Tree object

    df: bool, default True. This parameter is used when calling the function recursively. User should user df=True. 
        If True, return the distance matrix as a pandas Dataframe of shape (n_leaves, n_leaves). 
        If False, return two lists: a list of tuple [(leaf1, leaf2, distance), ...] and a dictionary {leaf: distance to root, ...}. 


    Returns:
    --------
    (Default) If df=True, return a pandas Dataframe of shape (n_leaves, n_leaves).

    If df=False, return two lists: a list of tuple [(leaf1, leaf2, distance), ...] and a dictionary {leaf: distance to root, ...}.
    """

    children = tree.children
    children_dist = [c.dist for c in children]

    all_pairwise = []
    all_root_dist = {}

    _root_dists = []

    for c in children:
        if c.is_leaf():
            pairwise = []
            root_dist = {c.name: 0}
        else:
            pairwise, root_dist = cal_dist_matrix_efficient(c, df=False)

        all_pairwise.extend(pairwise)
        _root_dists.append(root_dist)

        for k, v in root_dist.items():
            all_root_dist[k] = v+c.dist

    for (c1, c1_dist, c1_root_dist), (c2, c2_dist, c2_root_dist) in itertools.combinations(zip(children, children_dist, _root_dists), 2):
        for c1_leaf, c1_leaf_dist in c1_root_dist.items():
            for c2_leaf, c2_leaf_dist in c2_root_dist.items():
                all_pairwise.append(
                    [c1_leaf, c2_leaf, c1_leaf_dist+c2_leaf_dist+c1_dist+c2_dist])

    if df:
        df = {}
        for l, r, d in all_pairwise:
            try:
                df[l][r] = d
            except KeyError:
                df[l] = {r: d}
            try:
                df[r][l] = d
            except KeyError:
                df[r] = {l: d}

        df = pd.DataFrame(df)

        leaves = tree.get_leaf_names()
        df = df[leaves].loc[leaves, :]
        return df

    else:
        return all_pairwise, all_root_dist


class DataSplitter:
    """ Base clase for data spliters that split data into training sets and test sets.

    Child classes should implement at least two methods: __init__ and split. Implement generate_splits for efficiency.

    Methods
    -------
    generate_splits: should be the most used method. 

    split: generate a single split. Must be implemented by child classes. 

    """

    def __init__(self):
        raise NotImplementedError

    def split(self, samples):
        """ Generate one split. 

        Return: 
        -------
        A list of test set samples.
        """
        raise NotImplementedError

    def generate_splits(self, samples, n, **kwargs):  # TODO:  check for repeated splits
        """ Generate many splits by calling self.split repeatedly. Note that some splitters should implement this method for efficiency.

        Parameters:
        -----------
        samples: list of samples
        n: int. Number of splits to generate.

        Returns:
        --------
        A list of test set samples. 
        """
        return [self.split(samples, **kwargs) for _ in range(n)]


class LeaveOneOutSplitter(DataSplitter):
    def __init__(self):
        pass

    def generate_splits(self, samples, n=None, **kwargs):
        if n is None:
            n = len(samples)
        out = np.array([[s] for s in samples])
        if n < len(samples):
            out = np.random.choice(out, n, replace=False)
        return out


class RandomSplitter(DataSplitter):
    """ Random data split. 

    Parameters
    ----------
    test_set_ratio: float, default=0.2
        Ratio of test set size to the whole data set.

    """

    def __init__(self, test_set_ratio=0.2):
        self.test_set_ratio = test_set_ratio

    def split(self, samples):
        return np.random.choice(samples, size=int(len(samples)*self.test_set_ratio), replace=False)


class LargeTreeTraverseOOCSplitter(DataSplitter):
    """ Out-of-clade data split by iterating tree clades. Applicable to large trees. 
    
    Parameters
    ----------
    tree: ete3.Tree object

    test_set_range: tuple 
        Range of test set size ratio. 

    single_clades: list or None, default=None
        List of pre-computed single clades for the clade selection step. Precompute this to save time. 

    n_max_clade: int, default=2
        Max number of seperate clades in the test set. 

    prefer_small_clade: bool, default=False
        If true, single clade test sets are more likely to be selected. 
    
    growth_data: pd.Series 
        Binary growth data. Used to when min_zeros and min_ones are larger than 1. 
    
    min_zeros, min_ones: int, default=0
        Minimum sample numbers of zeros/ones in the test set. 
    
    time_out_iter: int or None, default=None
        Maximum iterations when trying to construct a test set.  

    """

    def __init__(self, tree,  
                 test_set_range=(0.2, 0.3),
                 single_clades=None,
                 n_max_clade=2, prefer_small_clade=False,
                 growth_data=None,min_zeros=0, min_ones=0, time_out_iter=None):

        self.tree = tree
        self.single_clades = single_clades
        self.growth_data = growth_data
        self.test_set_range = test_set_range
        self.min_zeros = min_zeros
        self.min_ones = min_ones
        
        if (bool(min_ones) or bool(min_zeros)) and growth_data is None:
            raise ValueError("Growth data is required if min_zeros or min_ones is larger than 1. ")

        self.prefer_small_clade = prefer_small_clade
        self.n_max_clade = n_max_clade
        if time_out_iter is None:
            self.time_out_iter = len(tree.get_leaves())*3
        else:
            self.time_out_iter = time_out_iter

    def split(self, samples, **kwargs):
        if self.single_clades is None:
            self.single_clades = self.compute_single_clades(self.tree, samples)
        i_iter = 0
        while True:
            if self.prefer_small_clade:
                clades = random.sample(
                    self.single_clades[1:], np.random.randint(self.n_max_clade)+1)
            else:
                clades = random.sample(self.single_clades, self.n_max_clade)
            test_samples = reduce(np.union1d, clades)
            if self.is_good_split(test_samples, samples):
                break
            i_iter += 1
            if self.time_out_iter and i_iter > self.time_out_iter:
                raise ValueError('Timeout in finding a good split.')
        return test_samples

    def compute_single_clades(self, tree, samples):
        clade_size_range = (int(
            self.test_set_range[0]*len(samples)), int(self.test_set_range[1]*len(samples)))
        
        tree=self.tree.copy()
        tree.prune(samples,preserve_branch_length=True)
        # tree = tree.copy()
        # tree_cleaned = False
        # while not tree_cleaned:
        #     tree_cleaned = True
        #     for l in tree.get_leaves():
        #         if l.name == '' or l.name not in samples:
        #             l.delete()
        #             tree_cleaned = False

        single_clades = [[]]  # including an empty clade
        for t in tree.traverse():
            if t.is_leaf():
                continue
            else:
                leaves = [l.name for l in t.get_leaves()]
                if len(leaves) > clade_size_range[1]:
                    continue
                else:
                    single_clades.append(leaves)

        return single_clades

    def is_good_split(self, test_samples, samples):
        if (len(test_samples)/len(samples)) < self.test_set_range[0] or (len(test_samples)/len(samples)) > self.test_set_range[1]:
            return False
        if self.min_zeros is not None and (self.growth_data[test_samples] == 0).sum() < self.min_zeros:
            return False
        if self.min_ones is not None and (self.growth_data[test_samples] == 1).sum() < self.min_ones:
            return False
        return True

# ============================================================
# Pipelines for running models
# ============================================================


def _pickleable_run_model(model, X_train, X_test, y_train, y_test, params, save_model):
    """ For multi-threading use in the PredictionPipeline class."""
    model.fit(X_train, y_train)
    score = model.score(X_test, y_test)
    params.update(model.get_params())
    params.update(score)
    params.update({'train_sample_size': len(
        params['train_samples']), 'test_sample_size': len(params['test_samples'])})
    if save_model:
        try: # If the model has a Pool attribute, remove it before pickling. 
            if model._p is not None:
                model._p = None
        except AttributeError:
            pass
        params['_model'] = model
    return params


class PredictionPipeline:
    """ The pipeline to evalutate prediction models. The pipeline generate training/test sets from a splitter, train/test models, and record prediction results. This class handles multi-threading and caching to save memories. 

    Typical usage:
    >>> pipe=PredictionPipeline(....) # Specifiy parameters
    >>> pipe.generate_splits(X,y) # Generate split samples
    >>> pipe.run() # The actual run

    Parameters:
    -----------
    Model : a BasePredictionModel class

    model_parames : dict. Default: None.
        Parameters for the model. 
        If None, use the default parameters.
        If not carbon-specific, the same parameters are used for all carbons.
        If carbon specific, format the dictionary as {c1: params_1, ...}

    threads : int. Default: None.
        Number of threads to use.  
        If >1, use the passed Pool object (p) or a new Pool object is created (if p == None). Passing a Pool object is recommented because opening multiple pools hurts performance significantly. 
            - If a Pool object is passed, this parameters is actually not used. Just put in a number larger than 1 to enable multi-threading.

    split_method : str, a DataSplitter class, or a dict. Default: 'random'.
        'random': random split (RandomSplitter). 
        'ooc': out-of-clade split (LargeTreeTraverseSplitter).
        (DEPRECATED) 'ooc_dist_mat': out-of-clade split using the distance matrix (DistMatOOCSplitter).
        'leave_one_out': leave-one-out split (LeaveOneOutSplitter).
        A DataSplitter class
        or, a dictionary: user-specified test set samples, specific to each carbon: {c1: [[s1,s2,...], [s1,s2,...], ...], ...}
    
    n_splits : int. Number of splits. Default: 10. 
        Passed into the Datasplitter.generate_splits method. Ignored if split_method is a dictionary.
    
    splitter_params : dict. Default: None.
        Parameters for the DataSplitter. 
        If None, use the default parameters for all carbons.
        If not carbon-specific, the same parameters are used for all carbons.
        If carbon specific, format the dictionary as {c1: params_1, ...}. 
    
    tree : ete3.Tree. Default: None.
        Required for out-of-clade splits. 
    
    carbons : list. Default: CARBONS.
        a list of carbons to evalutate models on. 
    
    save_models : bool. Default: False.
        Whether to save the trained models objects. 
        Saving model objects takes lots of memory/disk space. Discouraged for large datasets. 
    
    multithreading_batch : int. Number of split-train-test cycles passed to the multi-threading pool. If caching is enabled, this is also the caching batch size. Default: None. This is necessary for large datasets to save memory.
        If None, use 5 * threads.
    
    pass_pool_to_model : bool. Whether to pass the Pool object to the model. Used only for multi-threaded models. Default: False.
    
    p : multiprocessing.Pool. A Pool object to use for multi-threading. Default: None.
        If None, a new Pool object is created.
        Passing a pre-created Pool object is strongly recommended because opening multiple Pool objects hurts performance significantly.
    
    ff_results : str. Caching file path. Default: None. Caching is necessary for large datasets to save memory.
        If None, caching is disabled. 
        If not None, the pipeline checks if the files exists and if so, rename it to a non-existing file name. In the run, the pipeline pickles a batch of split-train-test cycle results (1 for single-threading, multithreading_batch for multi-threading) to the file. At the end of the run, the pipeline concatenates results as a DataFrame and saves it to the file.
    
    allow_empty_set : bool. Default: False.
        Whether to allow empty train or test sets. This is to fit FBAClassifier into the pipeline where this parameter is set True. For other models, this should be False. 
    """

    def __init__(self, Model, model_params=None, threads=1,
                 split_method='random', n_splits=10, splitter_params=None,
                 tree=None,
                 carbons=CARBONS,
                 save_models=False,
                 multithreading_batch=None,
                 #pass_pool_to_model=False,
                 p=None, 
                 ff_results=None,
                 allow_empty_set=False):  # verbose=True,
        self.carbons = carbons
        self.Model = Model
        if (self.Model == GreedyFeatureSelection) and threads > 1:
            print("Warning: do not use pipeline multi-threading with meta learning models. Use multi-threading within the model instead.")
            threads = 1
        if model_params is None:
            self.model_params = {c: {} for c in carbons}
        elif self.carbons[0] not in model_params:
            self.model_params = {c: model_params for c in carbons}
        else:
            self.model_params = model_params

        self.split_method = split_method
        self.n_splits = n_splits
        if splitter_params is None:
            self.splitter_params = {c: {} for c in carbons}
        elif self.carbons[0] not in splitter_params:
            self.splitter_params = {c: splitter_params for c in carbons}
        self.splitter_params = splitter_params
        self.tree = tree

        self.threads = threads
        if multithreading_batch is None:
            multithreading_batch = threads * 5
        self.multithreading_batch = multithreading_batch
        # self.pass_pool_to_model = pass_pool_to_model
        # if p is None:
        #     self._external_pool = False
        # else:
        #     self._external_pool = True
        if self.threads>1 and p is None:
            raise ValueError("A Pool object is required for multi-threading.")
        self._p = p

        # Caching
        if ff_results is not None and os.path.exists(ff_results):
            print(f"Warning: output file {ff_results} exists. ")
            i = 0
            _spl = ff_results.split('.')
            while os.path.exists(ff_results):
                ff_results = '.'.join(_spl[:-1])+'_'+str(i)+'.'+_spl[-1]
                i += 1
            print(f"Changed to {ff_results}. ")
        self.ff_results = ff_results
        self._results = []

        self.carbons = carbons
        self.save_models = save_models
        # self.verbose=verbose
        self.allow_empty_set = allow_empty_set

        # a list of sample labels for data splits. [(c, train_samples, test_samples),...]
        self.splits = None
        self.X = None
        self.y = None
        # a batch of data splits that is feeded into _pickleable_run_model for multi-threading. This is done in batch for save memory. [(c, X_train, X_test, y_train, y_test, params),...]
        self._data_split_batch = None

    def generate_splits(self, X, y):
        """ Generate splits (a dictinary of test set samples for each carbon) for the pipeline. Must call before run(). 

        Parameters
        ----------
        X : pandas.DataFrame. 
            The feature matrix of shape (n_samples, n_features). The index is sample labels and the columns are feature names.
        
        y : pandas.DataFrame. 
            The binary trait matrix of shape (n_samples, n_carbons). The index is sample labels and the columns are carbon names. Samples with NaN values are ignored for the specific carbons. 
        
        """
        self.X = X
        self.y = y
        if isinstance(self.split_method, type) or self.split_method in ['random', 'ooc', 'ooc_distmat', 'leave_one_out']:
            split_dict = {}
            for c in self.carbons:
                sample_mask = y[y[c].notna()].index.values
                if isinstance(self.split_method, type):
                    splitter = self.split_method(**self.splitter_params)
                elif self.split_method == 'random':
                    splitter = RandomSplitter(**self.splitter_params)
                elif self.split_method == 'ooc':
                    splitter = LargeTreeTraverseOOCSplitter(
                        tree=self.tree, growth_data=y[c][sample_mask], **self.splitter_params)
                elif self.split_method == 'leave_one_out':
                    splitter = LeaveOneOutSplitter(**self.splitter_params)
                clades = splitter.generate_splits(sample_mask, self.n_splits)
                split_dict[c] = clades
            self.splits = self._generate_splits_from_dict(X, y, split_dict)

        elif isinstance(self.split_method, dict):
            self.splits = self._generate_splits_from_dict(
                X, y, self.split_method)
        else:
            raise ValueError(f"Unknown split method.")

    def _generate_splits_from_dict(self, X, y, split_dict):
        """ Generate splits (a dictinary of test set samples for each carbon) for the pipeline. This method eliminates samples in split_dict with NaN y values."""
        splits = []
        # if self.verbose:
        #     _iter=tqdm(self.carbons,desc="Generating splits...")
        # else:
        #     _iter=self.carbons
        _iter = tqdm(self.carbons, desc="Generating splits...")
        for c in _iter:
            # filter out samples with no growth data
            sample_mask = y[y[c].notna()].index.values
            for test_samples in split_dict[c]:
                test_samples = np.intersect1d(test_samples, sample_mask)
                train_samples = np.setdiff1d(sample_mask, test_samples)
                if (len(test_samples) > 0 and len(train_samples) > 0) or (self.allow_empty_set):
                    splits.append((c, train_samples, test_samples))

        return splits

    def run(self):
        """ Run the pipeline. """
        if self.splits is None:
            raise ValueError("Please generate splits first.")
        _data_split_batch = []

        # if self.threads > 1 and self._p is None:
        #     self._p = Pool(self.threads)
        # elif self.pass_pool_to_model:
        #     model_params = self.model_params[self.carbons[0]]
        #     if 'threads' in model_params and model_params['threads'] > 1:
        #         self._p = Pool(model_params['threads'])
        #         for c in self.carbons:
        #             self.model_params[c]['p'] = self._p
        for i_split, (c, train_samples, test_samples) in enumerate(tqdm(self.splits, desc='Training models...')):
            X_train, X_test = self.X.loc[train_samples,
                                         :], self.X.loc[test_samples, :]
            y_train, y_test = self.y.loc[train_samples,
                                         c], self.y.loc[test_samples, c]
            if not self.allow_empty_set:
                if X_train.shape[0] == 0 or X_test.shape[0] == 0:
                    continue
                if len(y_train) == 0 or len(y_test) == 0:
                    continue
            param = [self.Model(**(self.model_params[c])),
                     X_train, X_test, y_train, y_test]
            param.append(
                {'carbon_name': c, 'train_samples': train_samples, 'test_samples': test_samples})
            param.append(self.save_models)
            if self.threads <= 1:  # single thread
                self._save_new_results([_pickleable_run_model(*param)])
            else:  # multi-thread
                _data_split_batch.append(param)
                if (len(_data_split_batch) >= self.multithreading_batch) or i_split == len(self.splits)-1:
                    result = self._p.starmap(
                        _pickleable_run_model, _data_split_batch)
                    self._save_new_results(result)
                    _data_split_batch = []
                else:
                    continue

        # if self._p is not None and (not self._external_pool):
        #     self._p.close()
        #     self._p = None
        return self._load_final_results()

    def _save_new_results(self, results):
        """ Save a batch of results to self._results (not caching) or disk (caching)."""
        if not isinstance(results, list):
            results = [results]

        if self.ff_results is None:
            self._results.extend(results)
        else:
            with open(self.ff_results, 'ab') as f:
                pickle.dump(results, f)

    def _load_final_results(self):
        """ Load final results from self._results (not caching) or disk (caching; also convert results to a DataFrame)."""
        if self.ff_results is None:
            results = self._results
        else:
            results = []
            with open(self.ff_results, 'rb') as f:
                while True:
                    try:
                        results.extend(pickle.load(f))
                    except EOFError:
                        break
        return pd.DataFrame(results)


def run_multiple_models(models, datasets, DIR_output, 
                        run_models=True, p=None,DIR_cache=None,DIR_results=None, 
                        concatenate_output=True,ff_results_all=None,
                        **kwargs):
    """ Run the split-train-test pipelines for multiple models on multiple datasets.

    Parameters
    ----------


    """
    if DIR_cache is None:
        DIR_cache=os.path.join(DIR_output,'cache')
        try:
            os.makedirs(DIR_cache)
        except FileExistsError:
            pass
    if DIR_results is None:
        DIR_results=os.path.join(DIR_output,'results')
        try:
            os.makedirs(DIR_results)
        except FileExistsError:
            pass

    if run_models:
        print("Running models...")
        for model_name, (Model, pipe_params) in models.items():
            for dataset_name,dataset in datasets.items():
                try:
                    print(f"Running {dataset_name} {model_name}")
                    ff_cache=os.path.join(DIR_cache,f'{dataset_name}_{model_name}.pk')
                    ff_results=os.path.join(DIR_results,f'{dataset_name}_{model_name}.pk')
                    if os.path.exists(ff_results):
                        print("Already exists. Skipping. ")
                        continue
                    # pamameter priority: pipe_params > kwargs > default
                    params={'tree':dataset['tree'], 'carbons': dataset['carbons'],'p':p, 'ff_results':ff_cache}
                    params.update(kwargs)
                    params.update(pipe_params)

                    pipe=PredictionPipeline(Model, **params)
                    pipe.generate_splits(dataset['ko_data'],dataset['growth_data'])
                    results=pipe.run()
                    with open(ff_results, 'wb') as f:
                        pickle.dump(results, f)
                    print(f"Finished {dataset_name} {model_name} ")
                except Exception as e: 
                    print(f"Failed to run {dataset_name} {model_name}")
                    traceback.print_exc()
    
        print("Running models done.")

    # Concatenate data
    if concatenate_output:
        print("Concatenating data...")
        if ff_results_all is None:
            ff_results_all=os.path.join(DIR_output,'results_all.pk')

        results_all=[]
        for model_name, (Model, pipe_params) in models.items():
            for dataset_name,dataset in datasets.items():
                try:
                    ff_results=os.path.join(DIR_results,f'{dataset_name}_{model_name}.pk')
                    results=pd.read_pickle(ff_results)
                    results['model']=model_name
                    results['dataset_name']=dataset_name
                    results_all.append(results)
                except Exception as e:
                    print(e)
                    print(f"Failed to load {dataset_name} {model_name}")

        results_all=pd.concat(results_all,axis=0,ignore_index=True)
        results_all.to_pickle(ff_results_all)
        print("Concatenating data done.")
        return results_all

# ==============================================================================
# Calculate statistics
# ==============================================================================



def compare_models(df, model_pairs,
                   seperate_by='carbon_name', model_key='model', metric='accuracy',
                   p_threshold=0.05,force_positive_t=True,
                    multi_testing_correction=False, correct_by='carbons'):
    """ Compare models and calcualte statistics from prediction results. 
    
    Parameters
    ----------
    df : pandas.DataFrame. 
        PredictionPipeline output. 

    model_pairs : list of tuples.
        Each tuple contains the names of models to be compared and a function to calculate statistics. The number of model names should equal to the number of argumetns of the function.

    seperate_by : str or a list of str. 
        keys to serperate the data by. Passed to df.groupby(). Default: 'carbon_name'. 

    model_key : str. 
        The column name for model names. Default: 'model'.

    metric : str. 
        The column name for metrics. Default: 'accuracy'.\

    p_threshold: None, or a float. Default: 0.05. 
        If None, no p_threshold parameter is passed to the function.
        If a float, the function is passed a p_threshold parameter. 
        
    multi_testing_correction : bool. Default: True. 
        If True, p_threshold is divided by the number of tests. 

    Returns
    -------
    A DataFrame of statistics. 
    """
    if isinstance(seperate_by, list) and len(seperate_by) == 1:
        seperate_by = seperate_by[0]
    gb = df.groupby(seperate_by)
    results = []
    if multi_testing_correction:
        p_threshold = p_threshold/len(gb)/len(model_pairs) # Bonferroni correction
    for name, group in tqdm(list(gb)):
        for _ in model_pairs:
            models = _[:-1]
            f = _[-1]
            # if p_threshold is not None:
            #     if multi_testing_correction:
            #         p_threshold = p_threshold/(len(gb)) # Bonferroni correction
            #     f = lambda *args: _[-1](*args, p_threshold=p_threshold)

            if not isinstance(seperate_by, list) and not isinstance(seperate_by, tuple):
                seperate_by = [seperate_by]
            if not isinstance(name, tuple):
                name = (name,)

            _prefix = '_'.join(models)

            try:
                arrs = [group[group[model_key] == model]
                        [metric].values for model in models]
                out = f(*arrs)
                # if isinstance(out, dict):
                for k, v in out.items():
                    result = dict(zip(seperate_by, name))
                    result['stat'] = _prefix+'_'+k
                    result['value'] = v
                    results.append(result)
                # if 'p' in out:
                #     significant= (out['p'] < p_threshold)
                #     if 't' in out and force_positive_t and out['t'] < 0:
                #         significant = False
                #         result['p']=0.999
                #     result = dict(zip(seperate_by, name))
                #     result['stat'] = _prefix+'_'+'significant'
                #     result['value'] = significant
                #     results.append(result)
                # else:
                #     result = dict(zip(seperate_by, name))
                #     result['stat'] = _prefix+'_'+'statistic'
                #     result['value'] = out
                #     results.append(result)
            except Exception as e:
                print(e)
    results=pd.DataFrame(results)
    # Mutiple testing correction
    def get_correct_pvalues(ps):
        return multipletests(ps, method=multi_testing_correction)[1]
    if multi_testing_correction:
        print("Correcting p-values...")
        if correct_by=="all":
            _ind=results[results['stat'].str.endswith('_p')].index.values
            results.loc[_ind, 'value'] = get_correct_pvalues(results.loc[_ind, 'value'].values)
        elif correct_by=="carbons":
            if (isinstance(seperate_by, list) and len(seperate_by)>1) or (isinstance(seperate_by,str)):
                raise ValueError
            seperate_by=seperate_by[0]
            for carbon_name in df[seperate_by].unique():
                _ind=results[(results['stat'].str.endswith('_p')) & (results[seperate_by]==carbon_name)].index.values
                results.loc[_ind, 'value'] = get_correct_pvalues(results.loc[_ind, 'value'].values)
        else:
            raise ValueError

    return results.sort_values(by='stat').pivot(index=seperate_by, columns='stat', values='value').reset_index()


def ttest(arr1, arr2):
    """ Calculate t-test statistics. 
    
    Parameters
    ----------
    arr1, arr2: array-like of length >1.
        Assume no NaN values.

    Returns
    -------
    A dict of t, p, and significant.
    """
    if len(arr1) == 1 or len(arr2) == 1:
        raise ValueError("array length must be greater than 1")
    t, p = ttest_ind(arr1, arr2, equal_var=False)
    return {'t': t, 'p': p}

def ttest_permutation(arr1, arr2, n_permutations=100000):
    """ Calculate t-test statistics and p-values from permutation. 
    
    Parameters
    ----------
    arr1, arr2: array-like of length >1. 
        Assume no NaN values.
  
    Returns
    -------
    A dict of t, p, and significant.
    """
    if len(arr1) == 1 or len(arr2) == 1:
        raise ValueError("array length must be greater than 1")
    t, p = ttest_ind(arr1, arr2, permutations=n_permutations, equal_var=False)
    return {'t': t, 'p': p}


def single_model_summary(arr):
    """ Calculate summary statistics for results of a single models.

    Parameters
    ----------
    arr: array-like. 
        Assume no NaN values.

    Returns
    -------
    A dict of mean, std, min, and max.
    """

    if len(arr) == 1:
        raise ValueError("array length must be greater than 1")
    return {'mean': arr.mean(), 'std': arr.std(), 'min': arr.min(), 'max': arr.max()}


def one_sample_test(arr1, arr2, sign='>'):
    """ When one array has one sample and the other array (a null model) has multiple samples. p-values is the proportion of samples in the null model that are greater/smaller (decided by the sign parameter) than the sample in the first array.

    Parameters
    ----------
    arr1: array-like of length 1. 

    arr2: array-like of length >1. 
        Assume no NaN values.

    sign: {'>','>=','<','<='}, default='>' 
        The p-value equals to mean(arr2 [sign] arr1).

    Returns
    -------
    A dict of p and significant.
    """
    if len(arr1) != 1 or len(arr2) == 1:
        raise ValueError(
            "First array length must be 1 and second array length must be greater than 1")

    if sign == '>':
        p = len(arr2[arr2 > arr1])/len(arr2)
    elif sign == '>=':
        p = len(arr2[arr2 >= arr1])/len(arr2)
    elif sign == '<':
        p = len(arr2[arr2 < arr1])/len(arr2)
    elif sign == '<=':
        p = len(arr2[arr2 <= arr1])/len(arr2)
    else:
        raise ValueError("sign must be one of ['>','>=','<','<=']")
    # if p==0:
    #     p=1/len(arr2) # avoid p=0
    return {'p': p}


# ============================================================
# Plotting functions
# ============================================================


def plot_tree_matrix(tree, matrix, cbar=True,
                     figsize=(15, 20), width_ratio=None,
                     titles=None,
                     yticks=None, vmin=0, vmax=1,
                     tree_params={}):
    """ Plot a tree next to matrices. Because ete3 does not support matplotlib, I used Bio.Phylo for a simple tree plot. 
    
    Parameters
    ----------
    tree: ete3.Tree. 
        The tree leaves and the matrix index must be identical. 

    matrix: a DataFrame or a list DataFrame. 

    cbar: bool or a list of bool of size (1+len(matrix),). Default: True.
        Whether to plot colorbar. 

    figsize: tuple. Default: (15,20).

    width_ratio: None, or a list of int of shape (1+len(matrix),). Default: determined by matrix shapes.
        The width ratio of the tree and the matrices. 

    titles: None or a list of str of shape (len(matrix),). Default: ['']*len(matrix).
        The titles of the matrices. 

    yticks: None/bool or a list of bool of shape (len(matrix),). Default: None. 
        Whether to plot yticks. 

    vmin, vmax: float. Default: 0, 1.
        Passed into sns.heatmap(...). 

    tree_params: dict. 
        Plotting parameters for Bio.Phylo plotting. 
    """

    tree_bio = Phylo.read(StringIO(tree.write()), format='newick')
    tree_bio.rooted = True
    leaves = tree.get_leaf_names()

    if isinstance(matrix, pd.DataFrame):
        matrix = [matrix]
    if isinstance(cbar, bool):
        cbar = [cbar]*len(matrix)

    if width_ratio is None:
        width_ratio = [2] + [m.shape[1] for m in matrix]

    if titles is None:
        titles = ['']*len(matrix)

    if yticks is None:
        yticks = [True]+[False]*(len(matrix)-1)
    elif isinstance(yticks, bool):
        yticks = [yticks]*len(matrix)

    for m in matrix:
        if len(np.intersect1d(leaves, m.index.values)) != len(leaves):
            raise ValueError("Tree and matrix do not have the same samples.")

    fig, axes = plt.subplots(
        1, 1+len(matrix), figsize=figsize, gridspec_kw={'width_ratios': width_ratio})

    # show_branch_support=False,show_scale=False,show_branch_length=False,
    Phylo.draw(tree_bio, axes=axes[0], do_show=False, **tree_params)
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[0].set_xlabel('')
    axes[0].set_ylabel('')
    axes[0].set_frame_on(False)
    axes[0].set_title("Tree")

    for m, cb, ax, title, ytick in zip(matrix, cbar, axes[1:], titles, yticks):
        m = m.loc[leaves, :]
        sns.heatmap(m, ax=ax, cbar=cb, vmin=vmin, vmax=vmax)
        if not ytick:
            ax.set_yticks([])
        ax.set_title(title)
        ax.set_ylabel('')

    fig.tight_layout()
    return fig, axes


def cluster_matrix(matrix, axis=0, metric=None, **kwargs):
    """ Hierarchical cluster a matrix. 
    
    """
    if axis == 1:
        return cluster_matrix(matrix.T, axis=0, metric=metric, **kwargs).T

    if metric is None:
        def metric(x, y):
            x = (x > 0).astype(int)
            y = (y > 0).astype(int)
            # -np.logical_and(np.isnan(x),np.isnan(y)).sum())
            return 1-(x == y).sum()/len(x)
    Z = hierarchy.linkage(matrix, metric=metric, **kwargs)
    return matrix.iloc[hierarchy.leaves_list(Z), :]


def plot_prediction_pipeline_results(results, growth_data, tree, carbons=CARBONS, metric='accuracy', tree_width=10, matrix_width=10):
    """ Quick visualization of prediction pipeline results. Plot tree, data split, growth data and model predictions all next to each other.


    
    """
    # data split - accuracy matrix
    matrices = []
    titles = []
    cbar = []
    width_ratio = [tree_width]

    for c, group in results.groupby('carbon_name'):
        if c not in carbons:
            continue

        split_matrix = pd.DataFrame(np.nan, index=growth_data.index, columns=[
                                    c]+list(range(group.shape[0])))
        split_matrix.iloc[:, 0] = growth_data.loc[:, c] # The column is the growth data. 
        for i, (_, row) in enumerate(group.iterrows()):
            split_matrix.loc[row['test_samples'], i] = row[metric]
        split_matrix.iloc[:, 1:] = cluster_matrix(
            split_matrix.iloc[:, 1:], axis=1, optimal_ordering=True).values
        matrices.append(split_matrix)
        titles.append(c+'_'+metric)
        cbar.append(True)
        width_ratio.append(matrix_width)

    fig, axes = plot_tree_matrix(tree, matrices, titles=titles, cbar=cbar, figsize=(
        5+5*len(carbons), 20), width_ratio=width_ratio)
    return fig


def plot_prediction_accuracy(result, ax=None, figsize=(10, 8), metric='accuracy', ylim=(-0.1, 1.1), hue=None):
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.gca()
    sns.violinplot(data=result, x='carbon_name', y=metric, ax=ax, hue=hue)
    sns.stripplot(data=result, x='carbon_name', y=metric, ax=ax, hue=hue)

    ax.grid()
    ax.set_ylim(*ylim)
    # x axis rotate 90 degree
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    return ax


def plot_model_comparison(models,
                          metric='accuracy', p_pairs=None, carbons=CARBONS, pvalue_threshold=0.05,
                          axes=None, figsize=(15, 25),):
    
    model_labels, null_labels = [], []
    models_concat = []

    for k, model in models.items():
        model = model[['carbon_name', metric]].copy()
        model['model_label'] = k
        if k.startswith("null"):
            null_labels.append(k)
            model['is_null'] = True

        else:
            model_labels.append(k)
            model['is_null'] = False
        models_concat.append(model)
    models_concat = pd.concat(models_concat, axis=0, ignore_index=True)
    models_concat = models_concat[models_concat['carbon_name'].isin(carbons)]

    def cal_p_value(arr1, arr2):
        if len(arr1) > 1 and len(arr2) > 1:
            return ttest_ind(arr1, arr2).pvalue
        elif len(arr2) == 1:
            arr1, arr2 = arr2, arr1

        if len(arr2) == 1:
            return np.nan

        else:
            arr2 = np.array(arr2)
            return np.sum(arr2 > arr1)/len(arr2)

    if p_pairs is None:
        p_pairs = [None]
    pvalues = []
    for c in carbons:
        for p_pair in p_pairs:
            if p_pair is None:  # all model - null pairs
                for model_label in model_labels:
                    for null_label in null_labels:
                        pvalues.append(
                            {'carbon_name': c, 'model1': model_label, 'model2': null_label, 'pvalue': np.nan})
            else:
                pvalues.append(
                    {'carbon_name': c, 'model1': p_pair[0], 'model2': p_pair[1], 'pvalue': np.nan})

    model_gb = models_concat.groupby(['model_label', 'carbon_name'])
    for d in pvalues:
        m1, m2 = model_gb.get_group((d['model1'], d['carbon_name'])), model_gb.get_group(
            (d['model2'], d['carbon_name']))
        try:
            d['pvalue'] = cal_p_value(m1[metric].values, m2[metric].values)
        except:
            pass
    pvalues = pd.DataFrame(pvalues)
    pvalues['model_pair'] = pvalues['model1']+" ; "+pvalues['model2']
    pvalues['significant'] = pvalues['pvalue'] < pvalue_threshold

    if axes is None:
        fig, axes = plt.subplots(len(carbons), 2, figsize=figsize, gridspec_kw={
                                 'width_ratios': [1, 0.5]})
    if not isinstance(axes[0], list):
        axes = [axes]
    for (ax0, ax1), (c, models_group), (_, pvalues_group) in zip(axes, models_concat.groupby('carbon_name'), pvalues.groupby('carbon_name')):

        if models_group.shape[0] > 1:
            sns.violinplot(x='model_label', y=metric, data=models_group,
                           ax=ax0, hue='is_null', order=[*model_labels, *null_labels])
        sns.stripplot(x='model_label', y=metric, data=models_group,
                      ax=ax0, hue='is_null', order=[*model_labels, *null_labels])
        ax0.set_title(c)
        ax0.grid()
        ax0.set_ylim(-0.1, 1.1)
        # disable color bar
        ax0.legend().set_visible(False)

        sns.scatterplot(x='model_pair', y='pvalue',
                        data=pvalues_group, ax=ax1, hue='significant')
        ax1.grid()

        # rotate xtick labels and center
        for tick in ax1.get_xticklabels():
            tick.set_rotation(30)
            tick.set_horizontalalignment("right")
        # ax1.semilogy()
        ax1.set_ylim(-0.1, 1.1)

    plt.gcf().tight_layout()

    return plt.gcf(), axes



def plot_fancy_model_comparison(df,model_pairs, 
                                stats=None,  metric='accuracy', multi_testing_correction=False,
                                hue_order=None,colors=None, rotate_x=True, 
                                pair_annotation=False, 
                                single_annotation=True, null_models=None, y_single_annotation=1.05,
                                show_null_model=True, 
                                legend_col=5,
                                # legend_title=None, legend_names=None,
                                **kwargs):
    """ Fancy model comparison plots, including proper coloring and p-value annotation. 
    
    Parameters
    ----------
    df: pd.DataFrame
        Result DataFrame from PredictionPipeline.

    model_pairs : list of tuples
        Pairs of models to be compared.
    
    stats: pd.DataFrame
        Statistics calculated from compare_models. If None, compute using model pairs using ttest_permutation. Pre-compute is recommended because permutation tests are slow. 
    
    multi_testing_correction: bool
        If True, apply multiple testing correction to p-values.
    
    hue_order: list or None, default=None
        hue_order in sns.catplot(). Only models in hue_order are plotted.
    
    colors: dict or None
        Dictionary of colors for models. Null models are automatically assigned grey. If None, assigned automatically.
    
    kwargs: keyward arugments passed to sns.catplot

    """
    


    if hue_order is None:
        hue_order=np.unique(np.array(model_pairs).flatten())
    df=df[df['model'].isin(hue_order )]

    if null_models is None:
        null_models=[m for m in df['model'].unique() if 'null' in m]
        print(null_models)
    
    if not show_null_model:
        df=df[~(df['model'].isin(null_models))]
        hue_order=[m for m in hue_order if m not in null_models]


    if colors is None:
        # 10 colors
        _models=[m for m in hue_order if 'null' not in m]
        colors = {model:color for model,color in zip(_models,COLORS[:len(_models)])}
    
    for model in hue_order:
        if 'null' in model and model not in colors:
            if model.startswith('null'):
                colors[model]='grey'
            elif model.startswith('identity_null'):
                colors[model]='silver'
    

    _catplot_kwargs=dict(data=df, x='carbon_name',y=metric,
                hue='model', hue_order=hue_order, palette=colors,
                kind='violin',cut=0,linewidth=0.75,inner='point',dodge=True,legend_out=False)
    _catplot_kwargs.update(kwargs)

    sns.catplot(**_catplot_kwargs)
                # height=4, aspect=5,
    sns.stripplot(data=df, x='carbon_name',y=metric,
                hue='model', hue_order=hue_order,dodge=True,color='black', legend=False,
                size=1)
    sns.move_legend(plt.gca(),loc='upper center', bbox_to_anchor=(0.5, 1.3), ncol=legend_col)

    ax=plt.gca()
    ax.yaxis.grid(color='gray')
    # ax.grid()
    if rotate_x:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    
    # For single data points, draw a colored line by by wizardry 

    half_violin=0.4/len(hue_order)
    def get_x_coord(i_x, i_hue):
        # wizardry
        x=i_c + half_violin * (-len(hue_order)+2*i_hue+1)
        return x

    gb=df.groupby(['carbon_name','model'])
    epsilon=0.03 
    single_coords={}

    for i_c,carbon in enumerate(ax.get_xticklabels()):
        for i_hue, model in enumerate(hue_order):
            try:
                group=gb.get_group((carbon.get_text(),model))
            except KeyError:
                pass
            if group.shape[0]==1:
                x=get_x_coord(i_c, i_hue) # wizardry
                y=group[metric].values[0]
                try:
                    single_coords[model]['x'].append(x)
                    single_coords[model]['y'].append(y)
                except KeyError:
                    single_coords[model]={'x':[x], 'y':[y]}
                            
                # ax.plot([x-half_violin+epsilon,x+half_violin-epsilon],[y,y],color=colors[model],linewidth=2)
    for model in single_coords:
        ax.plot(single_coords[model]['x'],single_coords[model]['y'],'X',color=colors[model],markersize=10)
    # ax.plot(single_xs,single_ys,'X',color=colors['model'],markersize=3)

    # Calculate statistics
    if stats is None:
        raise ValueError
        # stats=compare_models(df, 
        #             model_pairs=[
        #                 (*model_pair, ttest_permutation) for model_pair in model_pairs
        #             ],
        #             seperate_by='carbon_name',
        #             model_key='model',
        #             metric=metric,
        #             p_threshold=0.05,
        #             multi_testing_correction=False)
    try:
        stats=stats.set_index('carbon_name')
    except KeyError:
        pass

    if single_annotation: 
        for i_c,carbon in enumerate(ax.get_xticklabels()):
            for i_hue, model in enumerate(hue_order):
                if model in null_models:
                    continue
                significant=True
                for null_m in null_models:
                    if (('ooc' in null_m) and ('ooc' not in model)) or (('ooc' not in null_m) and ('ooc' in model)):
                        continue
                    significant=significant and (stats.at[carbon.get_text(),f'{model}_{null_m}_p']<=0.05)
                    try:
                        significant=significant and (stats.at[carbon.get_text(),f'{model}_{null_m}_t']>0)
                    except KeyError:   
                        pass          
                if significant:
                    x=get_x_coord(i_c, i_hue)
                    ax.text(x, y_single_annotation,'*',fontsize=20, ha='center',va='bottom',color='black')


    
    # Pair annotation via statannotation
    if pair_annotation:
        pairs=[]
        p_values=[]
        for c in df['carbon_name'].unique():
            for m1,m2 in model_pairs:
                pairs.append(((c,m1),(c,m2)))
                p_values.append(stats.at[c,f'{m1}_{m2}_p'])

        annot=Annotator(ax, pairs, **_catplot_kwargs)
        if multi_testing_correction:
            raise ValueError # This should be done before hand
            # comparison_correction='Bonferroni'
        else:
            comparison_correction=None
        annot.configure(test=None, comparisons_correction=comparison_correction).set_pvalues(p_values).annotate()
    
    
    plt.ylim(bottom=0, top=1.2)
    ax.set_yticks(ax.get_yticks()[ax.get_yticks()<=1])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')

    return plt.gcf(), stats

def cal_feature_importance(results, kos_data, highlight=None, trim=100):
    df_fi = []
    if highlight is None or isinstance(highlight, list) or isinstance(highlight, np.ndarray) or (isinstance(highlight, dict) and results['carbon_name'].iloc[0] not in highlight):
        highlight={c: highlight for c in results['carbon_name'].unique()}
    for c, group in results.groupby('carbon_name'):
        arr_fis = np.array(list(group['feature_importances'].values))
        features = kos_data.columns
        df = pd.DataFrame({'features': features, 'carbon_name': [
                          c]*len(features), 'fi_mean': arr_fis.mean(axis=0), 'fi_std': arr_fis.std(axis=0)})
        hl=highlight[c]
        df['highlight'] = ''
        if isinstance(hl, list) or isinstance(hl, np.ndarray):
            df['highlight'] = df['features'].isin(highlight).astype(int)
        elif isinstance(hl, dict):
            for key, values in hl.items():
                df.loc[df[df['features'].isin(values)].index.values, 'highlight'] = key

        if trim and trim < df.shape[0]:
            df = df.sort_values('fi_mean', ascending=False).iloc[:trim]
        df_fi.append(df.reset_index())

    df_fi = pd.concat(df_fi, axis=0, ignore_index=True)
    return df_fi


# ============================================================
# Synthetic data model
# ============================================================

class BinaryLogicTree:
    def __init__(self):
        self.is_leaf = False
        self.left, self.right = None, None
        self.name = ''

    def __str__(self):
        """ Intepretable logic string"""
        if self.is_leaf:
            return str(self.name)
        left_str = str(self.left)
        if not self.left.is_leaf:
            left_str = f'({left_str})'
        right_str = str(self.right)
        if not self.right.is_leaf:
            right_str = f'({right_str})'

        return f'{left_str}{self.name}{right_str}'

    def to_newick(self):
        """ Newick format string"""
        if self.is_leaf:
            return self.name
        else:
            return f'({self.left.to_newick()},{self.right.to_newick()}):{self.name}'

    def random_populate(self, leaves, min_split_ratio=None, p_and=None):
        leaves = np.unique(leaves).astype(str)
        if len(leaves) == 1:
            self.is_leaf = True
            self.name = leaves[0]
            self._need_parathesis = False
            return self

        if min_split_ratio is None:
            min_n_leaves = 1
        elif min_split_ratio >= 0.5:
            raise ValueError
        else:
            min_n_leaves = int(len(leaves)*min_split_ratio)
        min_n_leaves = max(1, min_n_leaves)
        n_left = random.randint(min_n_leaves, len(leaves)-min_n_leaves)
        left_leaves = np.random.choice(leaves, n_left, replace=False)
        right_leaves = np.setdiff1d(leaves, left_leaves)

        left_tree = BinaryLogicTree()
        left_tree.random_populate(
            left_leaves, min_split_ratio=min_split_ratio, p_and=p_and)
        right_tree = BinaryLogicTree()
        right_tree.random_populate(
            right_leaves, min_split_ratio=min_split_ratio, p_and=p_and)

        self.left = left_tree
        self.right = right_tree
        self._is_leaf = False
        if p_and is None:
            p = None
        else:
            p = [p_and, 1-p_and]
        self.name = np.random.choice(['&', '|'], p=p)

    def calculate(self, value_map):
        if self.is_leaf:
            return value_map[self.name]
        else:
            if self.name == '&':
                return self.left.calculate(value_map) and self.right.calculate(value_map)
            elif self.name == '|':
                return self.left.calculate(value_map) or self.right.calculate(value_map)
            else:
                raise ValueError


# ============================================================
# Organize data
# ============================================================

def finalize_data(ko_data, growth_data, tree,remove_prefix=False,
                min_zeros=None,min_ones=None, min_growth_data_samples=None):    
    """ Keep samples with complete KO, growth, and tree data, and organize to a dictionary format. 

    Parameters
    ----------
    ko_data: pd.DataFrame (n_samples, n_features)
        KO matrix. Indices are samples and columns are KOs. 
    
    growth_data: pd.DataFrame (n_samples, n_carbons)
        Growth data. Indices are samples and columns are carbon name.
    
    tree: ete3.Tree
        Phylogenetic tree. 
    
    remove_prefix: bool, default=False
        If true, remove "_zeqian", "_matti", and "_bacdive" from leaf names in tree.
    
    """

    tree=tree.copy()
    if remove_prefix is not None:
        for leaf in tree.get_leaves():
            leaf.name=leaf.name.replace("matti_","").replace("zeqian_","").replace("bacdive_","")

    # Remove samples with incomplete data
    samples=np.intersect1d(ko_data.index,growth_data.index)
    samples=np.intersect1d(samples,[node.name for node in tree.get_leaves()])

    # Remove carbons with not enough samples
    growth_data=growth_data.loc[samples,:]
    if min_growth_data_samples is not None:
        cs=growth_data.columns[(growth_data.count(axis=0))>=min_growth_data_samples].values
        growth_data=growth_data.loc[:,cs]

    if min_zeros is not None:
        cs=growth_data.columns[((growth_data==0).sum(axis=0))>=min_zeros].values
        growth_data=growth_data.loc[:,cs]

    if min_ones is not None:
        cs=growth_data.columns[((growth_data==1).sum(axis=0))>=min_ones].values
        growth_data=growth_data.loc[:,cs]
    growth_data=growth_data.dropna(axis=0, how='all')
    samples=np.intersect1d(ko_data.index,growth_data.index)
    samples=np.intersect1d(samples,[node.name for node in tree.get_leaves()])

    print(f"{len(samples)} samples: ", samples)

    ko_data=ko_data.loc[samples]
    growth_data=growth_data.loc[samples]
    tree.prune(samples, preserve_branch_length=True)
    return {'ko_data':ko_data, 'growth_data':growth_data,'tree':tree,'samples': samples, 'carbons':growth_data.columns.values}
