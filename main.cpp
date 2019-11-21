#include "common.h"
#include <vector>

void simulateStep(std::vector<Node>& new_nodes,
                  const std::vector<Node>& nodes,
                  const stepParams){

    int length = stepParams.length;
    int width = stepParams.width;
    int depth = stepParams.depth;
    int deltax = stepParams.deltax;

    // case of second or later step
    for (int k = 0; k < depth; k++){
        for (int j = 0; j < width; j++){
            for (int i = 0; i < length; i++){
                int index = i + length * (j + width * k);
                Node node = nodes[index];
                Node *left = index % length > 0 ? &nodes[index-1] : NULL;
                Node *right = index % length < length - 1 ? &nodes[index+1] : NULL;
                Node *up = (index/length) % width > 0 ? &nodes[index-length] : NULL;
                Node *down = (index/length) % width < width - 1 ? &nodes[index+length] : NULL;
                Node *in = index/(length*width) > 0 ? &nodes[index-length*width] : NULL;
                Node *out = index/(length*width) < depth -1 ? &nodes[index+length*width] : NULL;
                Node next_node;

                float dudx = node.get_dudx(left, right);
                float dvdy = node.get_dudx(up,down);
                float dwdz = node.get_dwdz(in,out);
                float drhodx, drhody, drhodz;
                // update rhoair
                if (left && right){
                    drhodx = (right.rho_air - left.rho_air) / (2 * deltax);
                }
                else if (left){
                    drhodx = (node.rho_air - left.rho_air) / (deltax);
                }
                else if (right){
                    drhodx = (right.rho_air - node.rho_air) / (deltax);
                }
                else{
                    drhodx = 0.f;
                }
                if (up && down){
                    drhody = (up.rho_air - down.rho_air) / (2 * deltay);
                }
                else if (down){
                    drhody = (node.rho_air - down.rho_air) / (deltay);
                }
                else if (up){
                    drhody = (up.rho_air - node.rho_air) / (deltay);
                }
                else{
                    drhody = 0.f;
                }
                if (in && out){
                    drhodz = (in.rho_air - out.rho_air) / (2 * deltaz);
                }
                else if (out){
                    drhodz = (node.rho_air - out.rho_air) / (deltaz);
                }
                else if (in){
                    drhodz = (in.rho_air - node.rho_air) / (deltaz);
                }
                else{
                    drhodz = 0.f;
                }
                drhodt = -(rho_curr*(dudx+dvdy+dwdz)+drhodx*(node.u+node.v+node.w));
                next_node.rho_air = node.rho_air + drhodt;
                // update rhofuel
                if (left && right){
                    drhodx = (right.rho_fuel - left.rho_fuel) / (2 * deltax);
                }
                else if (left){
                    drhodx = (node.rho_fuel - left.rho_fuel) / (deltax);
                }
                else if (right){
                    drhodx = (right.rho_fuel - node.rho_fuel) / (deltax);
                }
                else{
                    drhodx = 0.f;
                }
                if (up && down){
                    drhody = (up.rho_fuel - down.rho_fuel) / (2 * deltay);
                }
                else if (down){
                    drhody = (node.rho_fuel - down.rho_fuel) / (deltay);
                }
                else if (up){
                    drhody = (up.rho_fuel - node.rho_fuel) / (deltay);
                }
                else{
                    drhody = 0.f;
                }
                if (in && out){
                    drhodz = (in.rho_fuel - out.rho_fuel) / (2 * deltaz);
                }
                else if (out){
                    drhodz = (node.rho_fuel - out.rho_fuel) / (deltaz);
                }
                else if (in){
                    drhodz = (in.rho_fuel - node.rho_fuel) / (deltaz);
                }
                else{
                    drhodz = 0.f;
                }
                drhodt = -(rho_curr*(dudx+dvdy+dwdz)+drhodx*(node.u+node.v+node.w));
                next_node.rho_fuel = node.rho_fuel + drhodt;
                // update rhoco2
                if (left && right){
                    drhodx = (right.rho_co2 - left.rho_co2) / (2 * deltax);
                }
                else if (left){
                    drhodx = (node.rho_co2 - left.rho_co2) / (deltax);
                }
                else if (right){
                    drhodx = (right.rho_co2 - node.rho_co2) / (deltax);
                }
                else{
                    drhodx = 0.f;
                }
                if (up && down){
                    drhody = (up.rho_co2 - down.rho_co2) / (2 * deltay);
                }
                else if (down){
                    drhody = (node.rho_co2 - down.rho_co2) / (deltay);
                }
                else if (up){
                    drhody = (up.rho_co2 - node.rho_co2) / (deltay);
                }
                else{
                    drhody = 0.f;
                }
                if (in && out){
                    drhodz = (in.rho_co2 - out.rho_co2) / (2 * deltaz);
                }
                else if (out){
                    drhodz = (node.rho_co2 - out.rho_co2) / (deltaz);
                }
                else if (in){
                    drhodz = (in.rho_co2 - node.rho_co2) / (deltaz);
                }
                else{
                    drhodz = 0.f;
                }
                drhodt = -(rho_curr*(dudx+dvdy+dwdz)+drhodx*(node.u+node.v+node.w));
                next_node.rho_co2 = node.rho_co2 + drhodt;
                // update rhonox
                if (left && right){
                    drhodx = (right.rho_nox - left.rho_nox) / (2 * deltax);
                }
                else if (left){
                    drhodx = (node.rho_nox - left.rho_nox) / (deltax);
                }
                else if (right){
                    drhodx = (right.rho_nox - node.rho_nox) / (deltax);
                }
                else{
                    drhodx = 0.f;
                }
                if (up && down){
                    drhody = (up.rho_nox - down.rho_nox) / (2 * deltay);
                }
                else if (down){
                    drhody = (node.rho_nox - down.rho_nox) / (deltay);
                }
                else if (up){
                    drhody = (up.rho_nox - node.rho_nox) / (deltay);
                }
                else{
                    drhody = 0.f;
                }
                if (in && out){
                    drhodz = (in.rho_nox - out.rho_nox) / (2 * deltaz);
                }
                else if (out){
                    drhodz = (node.rho_nox - out.rho_nox) / (deltaz);
                }
                else if (in){
                    drhodz = (in.rho_nox - node.rho_nox) / (deltaz);
                }
                else{
                    drhodz = 0.f;
                }
                drhodt = -(rho_curr*(dudx+dvdy+dwdz)+drhodx*(node.u+node.v+node.w));
                next_node.rho_nox = node.rho_nox + drhodt;
            }
        }
    }

    return;
}

int main(int argc, char *argv[]){
    StartupOptions options = parseOptions(argc, argv);
    return 0;
}