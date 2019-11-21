#include "common.h"
#include <vector>

#define first_deriv(varname,prev,curr,next,field,delta)\
                if (prev && next){ \
                    varname = (next->field - prev->field) / (2.f*delta); \
                } \
                else if (prev){ \
                    varname = (curr.field - prev->field) / (delta);\
                } \
                else if (next){ \
                    varname = (next->field - curr.field) / (delta);\
                } \
                else{ \
                    varneam = 0.f; \
                }

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
                Node *out = index/(length*width) > 0 ? &nodes[index-length*width] : NULL;
                Node *in = index/(length*width) < depth -1 ? &nodes[index+length*width] : NULL;
                Node next_node;

                float dudx = node.get_dudx(left, right);
                float dvdy = node.get_dudx(up,down);
                float dwdz = node.get_dwdz(in,out);
                float drhodx, drhody, drhodz;
                // update rhoair
                first_deriv(drhodx,left,node,right,rho_air,deltax);
                first_deriv(drhody,down,node,up,rho_air,deltay);
                first_deriv(drhodz,out,node,in,rho_air,deltaz);
                drhodt = -(rho_curr*(dudx+dvdy+dwdz)+drhodx*(node.u+node.v+node.w));
                next_node.rho_air = node.rho_air + drhodt;
                // update rhofuel
                first_deriv(drhodx,left,node,right,rho_fuel,deltax);
                first_deriv(drhody,down,node,up,rho_fuel,deltay);
                first_deriv(drhodz,out,node,in,rho_fuel,deltaz);
                drhodt = -(rho_curr*(dudx+dvdy+dwdz)+drhodx*(node.u+node.v+node.w));
                next_node.rho_fuel = node.rho_fuel + drhodt;
                // update rhoco2
                first_deriv(drhodx,left,node,right,rho_co2,deltax);
                first_deriv(drhody,down,node,up,rho_co2,deltay);
                first_deriv(drhodz,out,node,in,rho_co2,deltaz);
                drhodt = -(rho_curr*(dudx+dvdy+dwdz)+drhodx*(node.u+node.v+node.w));
                next_node.rho_co2 = node.rho_co2 + drhodt;
                // update rhonox
                first_deriv(drhodx,left,node,right,rho_nox,deltax);
                first_deriv(drhody,down,node,up,rho_nox,deltay);
                first_deriv(drhodz,out,node,in,rho_nox,deltaz);
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