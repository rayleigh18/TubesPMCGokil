#ifndef PROCEDURE_H
#define PROCEDURE_H
#include "configuration.h"
#include "stdlib.h"
// Initiation Tab. always use for initiate array of component
void initiateResTab(resistor_tab *tab){
    tab->Neff = 0;
    tab->array = (resistor_t*)malloc((tab->Neff) * sizeof(resistor_t));
}

void initiateInTab(inductor_tab *tab){
    tab->Neff = 0;
    tab->array = (inductor_t*)malloc((tab->Neff) * sizeof(inductor_t));
}

void initiateCapTab(capacitor_tab *tab){
    tab->Neff = 0;
    tab->array = (capacitor_t*)malloc((tab->Neff) * sizeof(capacitor_t));
}

void initiateVolTab(voltage_source_tab *tab){
    tab->Neff = 0;
    tab->array = (voltage_source_t*)malloc((tab->Neff) * sizeof(voltage_source_t));
}

void initiateCurTab(current_source_tab *tab){
    tab->Neff = 0;
    tab->array = (current_source_t*)malloc((tab->Neff) * sizeof(current_source_t));
}

// push component to array of component
void addResToTab(resistor_tab *tab , resistor_t res){
    tab->Neff += 1;
    tab->array = (resistor_t*)realloc((tab->array) , (tab->Neff)*sizeof(resistor_t));
    (tab->array)[(tab->Neff) + 1] = res;
}

void addInToTab(inductor_tab *tab , inductor_t in){
    tab->Neff += 1;
    tab->array = (inductor_t*)realloc((tab->array) , (tab->Neff)*sizeof(inductor_t));
    (tab->array)[(tab->Neff) + 1] = in;
}

void addCapToTab(capacitor_tab *tab , capacitor_t cap){
    tab->Neff += 1;
    tab->array = (capacitor_t*)realloc((tab->array) , (tab->Neff)*sizeof(capacitor_t));
    (tab->array)[(tab->Neff) + 1] = cap;
}

void addVolToTab(voltage_source_tab *tab , voltage_source_t vol){
    tab->Neff += 1;
    tab->array = (voltage_source_t*)realloc((tab->array) , (tab->Neff)*sizeof(voltage_source_t));
    (tab->array)[(tab->Neff) + 1] = vol;
}

void addCurToTab(current_source_tab *tab , current_source_t cur){
    tab->Neff += 1;
    tab->array = (current_source_t*)realloc((tab->array) , (tab->Neff)*sizeof(current_source_t));
    (tab->array)[(tab->Neff) + 1] = cur;
}

// make Component
resistor_t makeRes(double value, int node1 , int node2){
    resistor_t temp;
    temp.value = value;
    temp.node1 = node1;
    temp.node2 = node2;

    return temp;
}

inductor_t makeIn(double value, double last_cur, double node1, double node2){
    inductor_t temp;
    temp.value = value;
    temp.last_curr = last_cur;
    temp.node1 = node1;
    temp.node2 = node2;

    return temp;
}

capacitor_t makeCap(double value, double last_vol, double node1, double node2){
    capacitor_t temp;
    temp.value = value;
    temp.last_volt = last_vol;
    temp.node1 = node1;
    temp.node2 = node2;

    return temp;
}

current_source_t makeCurr(double value, int nodePos, int nodeNeg){
    current_source_t temp;
    temp.value = value;
    temp.nodePos = nodePos;
    temp.nodeNeg = nodeNeg;

    return temp;
}

voltage_source_t makeVolt(double value, int nodePos, int nodeNeg){
    voltage_source_t temp;
    temp.value = value;
    temp.nodePos = nodePos;
    temp.nodeNeg = nodeNeg;

    return temp;
}

#endif