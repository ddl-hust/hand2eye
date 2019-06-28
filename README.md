

We use a stationary camera (Camera mount out of robot)

- T_t_cal: transform from robot tool to calibration board (fixed)
- T_b_t: transform from robot base to tool
- T_c_cal: transform from camera to calibration board
- T_b_c: transform from robot base to camera (fixed). You need to calibrate this transform

Data is in the form of [x, y, z, ox, oy, oz, ow], position + quaternion

```
make Main  # build
./Main     # run
```

2019/6/24 11:09:05 
修改了hand2eye.cpp 轴角相加方式，不再理解为对应旋转矩阵乘法，只是简单的矢量相加。得到的标定矩阵我觉得ok,旋转矩阵的误差ok.
现在可以开始写说明文档了。
