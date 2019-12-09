import random
import math

filename = "input/uniform_1-1-10.txt";
deltax = 0.1;
deltay = 0.1;
deltaz = 0.01;
deltat = 0.000001;
numIterations = 1000;
length = 1;
width = 1;
depth = 10;
sparklen = 100;
sparkStart = 200;
sparkEnd = 300;

def rho_o2(x, y, z):
    return 0.2854;

def rho_n2(x, y, z):
    return 0.9396;

def rho_fuel(x, y, z):
    return 1.225 / 13.6;

def rho_co2(x, y, z):
    return 0.;

def rho_nox(x, y, z):
    return 0.;

def rho_h2o(x, y, z):
    return 0.;

def temperature(x, y, z):
    return 375.;

def u(x, y, z):
    lo = -.1 #if x < length/2 else -.2;
    hi = .1 #if x > length/2 else .2;
    return 0.# if x == 0 or x == length - 1 else random.uniform(lo,hi);

def v(x, y, z):
    lo = -.1 #if x < length/2 else -.2;
    hi = .1 #if x > length/2 else .2;
    return 0.# if y == 0 or y == width - 1 else random.uniform(lo,hi);

def w(x, y, z):
    # lo = -.15 + (-.85 - .15) * z / depth;
    # return random.uniform(lo,-0.1);
    return 0.# if z == 0 or z == depth - 1 else random.uniform(-.1,.1);
def spark(i):
    y = width // 2;
    x = length//4 + i * (length//2) / sparklen;
    z = -(4*depth/(3*sparklen**2)) * i**2 + (4 * depth/ (3 * sparklen)) * i;
    return x + y*length + z * length*width;

class Node(object):
    def __init__(self):
        pass

def calculate_pressure(node):
    nV = 0;
    nV += node.rho_o2 / 0.0319998;
    nV += node.rho_n2 / 0.0280134;
    nV += node.rho_fuel / 0.11426;
    nV += node.rho_co2 / 0.044098;
    nV += node.rho_nox / 0.0300061;
    nV += node.rho_h2o / 0.018015;
    return nV * 8.31 * node.temperature;
def temperature_to_internal_energy(temperature):
    return 268000 + (temperature - 375) /0.723;
def temperature_to_viscosity(temp):
    return 0.000001458 * math.sqrt(temp*temp*temp) / (temp + 110.4);
def temperature_to_conductivity(temp):
    return 1.3965 + (temp - 375.) * (-0.000056);

def main():
    f = open(filename,"w+");
    print("opened %s\n" % filename);
    f.write("%f\n" % deltax);
    f.write("%f\n" % deltay);
    f.write("%f\n" % deltaz);
    f.write("%f\n" % deltat);
    f.write("%d\n" % numIterations);
    f.write("%d\n" % length);
    f.write("%d\n" % width);
    f.write("%d\n" % depth);

    prev = 0.;

    for k in range(depth):
        for j in range(width):
            for i in range(length):
                node = Node();
                node.rho_o2 = rho_o2(i,j,k);
                node.rho_n2 = rho_n2(i,j,k);
                node.rho_fuel = rho_fuel(i,j,k);
                node.rho_co2 = rho_co2(i,j,k);
                node.rho_nox = rho_nox(i,j,k);
                node.rho_h2o = rho_h2o(i,j,k);
                node.temperature = temperature(i,j,k);
                node.pressure = calculate_pressure(node);
                node.internal_energy = temperature_to_internal_energy(node.temperature);
                node.viscosity = temperature_to_viscosity(node.temperature);
                node.conductivity = temperature_to_conductivity(node.temperature);
                node.u = u(i,j,k);
                node.v = v(i,j,k);
                node.w = w(i,j,k);
                f.write("%f " % node.rho_o2);
                f.write("%f " % node.rho_n2);
                f.write("%f " % node.rho_fuel);
                f.write("%f " % node.rho_co2);
                f.write("%f " % node.rho_nox);
                f.write("%f " % node.rho_h2o);
                f.write("%f " % node.pressure);
                f.write("%f " % node.temperature);
                f.write("%f " % node.viscosity);
                f.write("%f " % node.internal_energy);
                f.write("%f " % node.conductivity);
                f.write("%f " % node.u);
                f.write("%f " % node.v);
                f.write("%f\n" % node.w);
        percentage = 100 * k / depth;
        if (percentage - prev >= 10.):
            print("%f %% complete" % percentage)
            prev = percentage;
    f.write("%d\n" % sparklen);
    f.write("%d\n" % sparkStart);
    f.write("%d\n" % sparkEnd);
    for i in range(sparklen):
        f.write("%d\n" % spark(i));
    f.close();
    print("done\n");
    return 0;

main();