
import os
import subprocess
import urllib
import irlib
import numpy
import pytest

@pytest.mark.tryfirst()
def test_initialize_test_environment():
    """ Copy the test data locally """
    if not os.path.isdir("tests/data"):
        os.mkdir("tests/data")
    if not os.path.isfile("tests/data/test_argentiere.h5"):
        urllib.urlretrieve("https://dl.dropboxusercontent.com/u/375008/"
                           "irlib_test_data/test_argentiere.h5",
                           "tests/data/test_argentiere.h5")
    if not os.path.isfile("tests/data/test_milne.h5"):
        urllib.urlretrieve("https://dl.dropboxusercontent.com/u/375008/"
                           "irlib_test_data/test_milne.h5",
                           "tests/data/test_milne.h5")
    assert os.path.isfile("tests/data/test_argentiere.h5")
    assert os.path.isfile("tests/data/test_milne.h5")
    return

@pytest.fixture
def argenfnm():
    return "tests/data/test_argentiere.h5"

@pytest.fixture
def milnefnm():
    return "tests/data/test_milne.h5"

def test_h5_add_utm_argentiere(argenfnm):
    print(os.getcwd())
    print(os.listdir(os.getcwd()))
    print(os.listdir("tests/"))
    print(os.listdir("tests/data/"))
    stem, ext = os.path.splitext(argenfnm)
    utmfnm = stem + "_utm" + ext
    ret = subprocess.call(["h5_add_utm.py", argenfnm, utmfnm, "--swap_lon"])
    assert ret == 0
    S = irlib.Survey(utmfnm)
    L = S.ExtractLine(0)
    assert all(z==32 for z in L.metadata.zones)
    nan = numpy.nan
    n = numpy.array([ 346370.982676,  346370.982676,  346349.151958,  346349.15472 ,
                      346348.999665,  346348.999665,            nan,  346348.999665,
                      346348.41812 ,  346348.41812 ,            nan,  346345.003619,
                      346343.414045,  346343.5691  ,  346343.5691  ,  346343.5691  ,
                      346343.5691  ,  346343.5691  ,  346343.5691  ,  346343.5691  ,
                      346343.5691  ,  346343.5691  ,  346343.5691  ,  346343.5691  ,
                                nan,  346341.057314,  346341.057314,  346340.223679,
                      346340.223679,  346339.098926,  346339.098926,  346332.810278,
                      346332.810278,  346332.909735,  346332.909735,  346330.605529,
                      346324.740146,  346324.740146,  346323.521499,  346323.521499,
                      346318.935591,  346318.935591,  346316.124672,  346316.124672,
                                nan,            nan,  346308.237741,  346308.237741,
                      346307.511916,  346308.349951,  346308.349951,  346308.241722,
                      346308.241722,  346309.475826,  346309.475826,  346311.759659,
                      346311.759659,  346311.681825,  346311.681825,  346314.763159,
                      346316.892246,  346316.892246,  346316.96977 ,  346316.96977 ,
                      346316.96977 ,  346316.96977 ,            nan,            nan,
                      346316.96977 ,  346316.96977 ,  346316.892246,  346316.892246,
                      346316.892246,  346316.892246])
    nanidx = numpy.isnan(n)
    assert numpy.allclose(numpy.asarray(L.metadata.eastings)[~nanidx], n[~nanidx])
    return

def test_h5_add_utm_milne(milnefnm):
    stem, ext = os.path.splitext(milnefnm)
    utmfnm = stem + "_utm" + ext
    ret = subprocess.call(["h5_add_utm.py", milnefnm, utmfnm])
    assert ret == 0
    S = irlib.Survey(utmfnm)
    L = S.ExtractLine(1)

    nan = numpy.nan
    n = numpy.array([ 9161133.36736,  9161133.36736,  9161133.36736,  9161133.36745,
                      9161133.36761,  9161133.36761,  9161133.36745,  9161133.36736,
                      9161133.3672 ,  9161133.3672 ,  9161133.3672 ,  9161133.36704,
                      9161126.33162,  9161118.19976,            nan,  9161105.76789,
                      9161100.78348,  9161097.71067,  9161095.8658 ,            nan,
                      9161093.95717,  9161094.01193,  9161094.62972,  9161096.4692 ,
                      9161099.86058,  9161104.24245,  9161107.0551 ,  9161111.66542,
                      9161118.3973 ,  9161125.35879,  9161129.67392,  9161129.90568,
                      9161133.86133,  9161139.71231,  9161145.12097,  9161149.97233,
                      9161153.8066 ,  9161156.41596,  9161159.46658,  9161159.59457,
                      9161163.20475,  9161167.69076,  9161171.95208,  9161173.63734,
                      9161173.56501,  9161175.60739,            nan,  9161178.1333 ,
                      9161180.50754,  9161182.65744,  9161185.03814,  9161188.07447,
                      9161191.79024,  9161194.27968,  9161198.55899,            nan,
                      9161203.08675,  9161207.25111,  9161211.97951,            nan,
                      9161219.97287,  9161225.04166,  9161230.32451,  9161235.93926,
                      9161240.99503,  9161251.44201,  9161251.89183,  9161256.72228,
                      9161262.00746,  9161266.84627,  9161271.5694 ,            nan,
                      9161281.12757,  9161283.16876,  9161283.73181,  9161285.64676,
                      9161289.67984,  9161292.93198,  9161294.72435,  9161296.62431,
                                nan,  9161297.51707,  9161299.97351,  9161300.53191,
                      9161300.53191,            nan,  9161303.10884,  9161306.69693,
                      9161310.2838 ,  9161315.08688,  9161318.10779,  9161322.02093,
                      9161326.05301,  9161331.5344 ,  9161335.79309,  9161339.72034,
                      9161343.99263,  9161348.4838 ,  9161352.42586,  9161358.60804,
                      9161359.17132,  9161363.56273,  9161367.83803,  9161371.44716,
                      9161371.82944,  9161369.09333,            nan,  9161361.14841,
                      9161358.84448,  9161356.77308,  9161354.48075,  9161353.75103,
                                nan,  9161354.62899,  9161353.67833,  9161352.61849,
                      9161352.11253,  9161352.05144,  9161354.67132,            nan,
                      9161357.74864,  9161358.31425,  9161358.324  ,  9161359.26338,
                      9161360.07601,  9161365.70368,            nan,  9161375.94677,
                      9161379.22823,  9161383.50879,  9161387.7869 ,  9161391.61551,
                      9161394.44281])
    nanidx = numpy.isnan(n)
    assert numpy.allclose(numpy.asarray(L.metadata.northings)[~nanidx], n[~nanidx])
    assert all(z==17 for z in L.metadata.zones)

