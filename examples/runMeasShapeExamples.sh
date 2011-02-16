EX=$EXTEND_SHAPE_HSM_DIR/examples
DATA=$EXTEND_SHAPE_HSM_DIR/tests/data
$EX/meas_shape $DATA/image.0.fits $DATA/psf.0.fits 35.888 19.845 35.012 REGAUSS > myexample.regauss.out
$EX/meas_shape $DATA/image.2.fits $DATA/psf.2.fits 19.440 25.047 35.934 REGAUSS >> myexample.regauss.out
$EX/meas_shape $DATA/image.4.fits $DATA/psf.4.fits 8.740 11.920 35.155 REGAUSS >> myexample.regauss.out
$EX/meas_shape $DATA/image.6.fits $DATA/psf.6.fits 20.193 38.930 35.111 REGAUSS >> myexample.regauss.out
$EX/meas_shape $DATA/image.8.fits $DATA/psf.8.fits 57.940 27.73 35.165 REGAUSS >> myexample.regauss.out
$EX/meas_shape $DATA/image.0.fits $DATA/psf.0.fits 35.888 19.845 35.012 BJ > myexample.bj.out
$EX/meas_shape $DATA/image.2.fits $DATA/psf.2.fits 19.440 25.047 35.934 BJ >> myexample.bj.out
$EX/meas_shape $DATA/image.4.fits $DATA/psf.4.fits 8.740 11.920 35.155 BJ >> myexample.bj.out
$EX/meas_shape $DATA/image.6.fits $DATA/psf.6.fits 20.193 38.930 35.111 BJ >> myexample.bj.out
$EX/meas_shape $DATA/image.8.fits $DATA/psf.8.fits 57.940 27.73 35.165 BJ >> myexample.bj.out
$EX/meas_shape $DATA/image.0.fits $DATA/psf.0.fits 35.888 19.845 35.012 LINEAR > myexample.linear.out
$EX/meas_shape $DATA/image.2.fits $DATA/psf.2.fits 19.440 25.047 35.934 LINEAR >> myexample.linear.out
$EX/meas_shape $DATA/image.4.fits $DATA/psf.4.fits 8.740 11.920 35.155 LINEAR >> myexample.linear.out
$EX/meas_shape $DATA/image.6.fits $DATA/psf.6.fits 20.193 38.930 35.111 LINEAR >> myexample.linear.out
$EX/meas_shape $DATA/image.8.fits $DATA/psf.8.fits 57.940 27.73 35.165 LINEAR >> myexample.linear.out
$EX/meas_shape $DATA/image.0.fits $DATA/psf.0.fits 35.888 19.845 35.012 KSB > myexample.ksb.out
$EX/meas_shape $DATA/image.2.fits $DATA/psf.2.fits 19.440 25.047 35.934 KSB >> myexample.ksb.out
$EX/meas_shape $DATA/image.4.fits $DATA/psf.4.fits 8.740 11.920 35.155 KSB >> myexample.ksb.out
$EX/meas_shape $DATA/image.6.fits $DATA/psf.6.fits 20.193 38.930 35.111 KSB >> myexample.ksb.out
$EX/meas_shape $DATA/image.8.fits $DATA/psf.8.fits 57.940 27.73 35.165 KSB >> myexample.ksb.out
$EX/meas_shape $DATA/image.0.fits $DATA/psf.0.fits 35.888 19.845 35.012 SHAPELET8,8 > myexample.shapelet8.out
$EX/meas_shape $DATA/image.2.fits $DATA/psf.2.fits 19.440 25.047 35.934 SHAPELET8,8 >> myexample.shapelet8.out
$EX/meas_shape $DATA/image.4.fits $DATA/psf.4.fits 8.740 11.920 35.155 SHAPELET8,8 >> myexample.shapelet8.out
$EX/meas_shape $DATA/image.6.fits $DATA/psf.6.fits 20.193 38.930 35.111 SHAPELET8,8 >> myexample.shapelet8.out
$EX/meas_shape $DATA/image.8.fits $DATA/psf.8.fits 57.940 27.73 35.165 SHAPELET8,8 >> myexample.shapelet8.out
