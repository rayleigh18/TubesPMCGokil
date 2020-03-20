#ifndef INPUT_H
#define INPUT_H

#include "stdio.h"
#include "../Configuration/configuration.h"
resistor_t askResistor(){
    resistor_t temp;
    printf("Masukkan nilai Resistor :");
    scanf("%lf", &(temp.value));
    printf("Masukkan node :");
    scanf("%d", &(temp.node1));
    printf("Masukkan node :");
    scanf("%d", &(temp.node2));
    return temp;  
}

capacitor_t askCapacitor(){
    capacitor_t temp;
    printf("Masukkan nilai Capacitor :");
    scanf("%lf", &(temp.value));
    printf("Masukkan node :");
    scanf("%d", &(temp.node1));
    printf("Masukkan node :");
    scanf("%d", &(temp.node2));
    printf("Masukkan nilai Tegangan saat t = t0 :");
    scanf("%lf", &(temp.last_volt));
    return temp;  
}

inductor_t askInductor(){
    inductor_t temp;
    printf("Masukkan nilai Induktor :");
    scanf("%lf", &(temp.value));
    printf("Masukkan node :");
    scanf("%d", &(temp.node1));
    printf("Masukkan node :");
    scanf("%d", &(temp.node2));
    printf("Masukkan nilai Arus saat t = t0 :");
    scanf("%lf", &(temp.last_curr));
    return temp;  
}

voltage_source_t askVoltageSource(){
    voltage_source_t temp;
    printf("Masukkan nilai Tegangan :");
    scanf("%lf", &(temp.value));
    printf("Masukkan nilai Node Positif :");
    scanf("%d", &(temp.nodePos));
    printf("Masukkan nilai Node Negatif :");
    scanf("%d", &(temp.nodeNeg));
    
    return temp;
}

current_source_t askCurrentSource(){
    current_source_t temp;
    printf("Masukkan nilai Arus :");
    scanf("%lf", &(temp.value));
    printf("Masukkan nilai Node Positif :");
    scanf("%d", &(temp.nodePos));
    printf("Masukkan nilai Node Negatif :");
    scanf("%d", &(temp.nodeNeg));

    return temp;
}

void pushResToNodeArray(resistor_t res,int num_in_array ,node_tab *node_list){
    int i = 0;
    int found1 = 0;
    int found2 = 0;
    for (i = 0; i < node_list->Neff; i++){
        int tempName = (node_list->array)[i].name;
        if (tempName == res.node1){
            found1 = 1;
            addIntegerToTab(&((node_list->array)[i].res_list) , num_in_array);
        }
        else if (tempName == res.node1){
            found2 = 1;
            addIntegerToTab(&((node_list->array)[i].res_list) , num_in_array);
        }
    }

    if (found1 != 1){
        node_t tempNode = makeNode(res.node1);
    }
}

#endif