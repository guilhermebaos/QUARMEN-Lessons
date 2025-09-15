# Synchronization 

## Master
We can specify a region that will only be executed by the master thread:


### Use Cases
This command can be used for example to:
- Print output to the user.
- Use the values of the various threads to produce one final result.



# Work-Sharing

## For Loops
We can work with for loops in two ways.

### Inside an Existing Parallel Block
We can create a for loop inside a previously created parallel block:


## False Sharing
When we create a variable we store it physically on the memory, but different threads have privilleged access to certain areas of memory, hence when we share memory accross threads we might 
