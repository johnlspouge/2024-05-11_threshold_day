#!/usr/bin/env python
"""
Heap_Dictionary sorts a dictionary on its keys with a heap.
For performance, it can resolve_collision a key collision with an arbitrary operation.
"""

# Events stores epidemic events like infection and compartmental decay.
#   All events should pass through the object by push() and pop().
#   pop() and set_current() (called to initialize the output with t=0) trigger the callback.

import heapq

# Heap_Dictionary is a heap of dictionaries sorted by keys.
#   Heap_Dictionary can resolve key collisions with set_resolve_collision() for resolve_collision(k, self.dictionary, dictionary).
#   Heap_Dictionary maintains a current key, which may or may not be in the heap.
#   The current key can tracks the minimum element of the heap.
class Heap_Dictionary(): 
    def __init__(self): 
        self.keys = []
        self.dictionary = {}
    def clear(self):
        self.keys = []
        self.dictionary = {}
        return self
    # Returns the dictionary for read-only. Otherwise, consistency with the heap is lost.
    def get_dictionary(self):
        return self.dictionary 
    # Sets the resolution of collisions of the same key with set_resolve_collision(key, self.d0, d).
    def set_resolve_collision(self, resolve_collision=None):
        self.resolve_collision = resolve_collision 
        return self
    def is_empty(self):
        return not bool(self.keys)
    # Peeks at the next key. Returns None if the keys are empty
    def peek(self): 
        if self.keys:
            return self.keys[0]
        return None
    # Pushes a dictionary onto the heap, never with a callback.
    def push(self, dictionary): # period2transition = {time:transition}. 
        for k,v in dictionary.items():
            if k in self.dictionary: # The key is already on the heap. Only the dictionary requires attention.
                if hasattr(self, "resolve_collision"): 
                    self.resolve_collision(k, self.dictionary, dictionary)
                else:
                    self.dictionary[k] = v
            else:  # new key
                heapq.heappush(self.keys, k)
                self.dictionary[k] = v
        assert len(self.keys) == len(self.dictionary)
        return self
    # Pops the next dictionary value.
    def pop(self):
        if self.is_empty():
            return None
        key = heapq.heappop(self.keys)
        value = self.dictionary[key]
        del self.dictionary[key]
        assert len(self.keys) == len(self.dictionary)
        return value
    # Returns a sorted list of keys. The log(len(self.keys())) operation is for debugging only.
    def get_sorted_list_of_keys(self):
        return sorted(self.dictionary.keys())

def test_Heap_Dictionary():
    hd = Heap_Dictionary()  
    def resolve_collision(k, d0, d):
        d0[k] += d[k]
    hd.set_resolve_collision(resolve_collision)
    hd.push({1: 4})
    hd.push({0: 1}).push({1:-2}) # 4-2 == 2, as from resolve_collision
    hd.push({2: 3})
    assert hd.dictionary == {k:(k+1) for k in range(0,3)}
    for i,d0 in enumerate([{0:1,1:2,2:3},{1:2,2:3},{2:3}]):
        assert hd.dictionary == d0
        assert hd.pop() == i + 1 # Pops the dictionary item with the smallest key.
    assert hd.is_empty()
    
def main():
    test_Heap_Dictionary()
    
if __name__ == "__main__":
    main()
