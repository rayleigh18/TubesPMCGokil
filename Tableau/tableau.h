void getInstance(
  component *component_array,
  int component_array_length,
  int ground,
  int *node_array,
  int node_array_length,
  double **instance,
  int *instance_length
);

void printInstance(
  double *instance,
  int instance_length,
  component *component_array,
  int component_array_length,
  int *node_array,
  int node_array_length,
  int ground
);