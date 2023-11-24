import time

# O(1) - Constant Time
def is_first_element_null(elements):
    """Check if the first element of the list is None."""
    return elements[0] is None

# O(n) - Linear Time
def find_max(elements):
    """Find the maximum element in a list."""
    max_value = elements[0]
    for num in elements:
        if num > max_value:
            max_value = num
    return max_value

# O(n^2) - Quadratic Time
def bubble_sort(elements):
    """Sort a list using the bubble sort algorithm."""
    n = len(elements)
    for i in range(n):
        for j in range(0, n-i-1):
            if elements[j] > elements[j+1]:
                elements[j], elements[j+1] = elements[j+1], elements[j]
    return elements

# O(2^n) - Exponential Time
def fibonacci(n):
    """Calculate the nth Fibonacci number."""
    if n <= 1:
        return n
    else:
        return fibonacci(n-1) + fibonacci(n-2)

# O(n!) - Factorial Time
def permute(elements):
    """Generate all permutations of a list."""
    if len(elements) == 1:
        return [elements]
    permutations = []
    for i in range(len(elements)):
        for p in permute(elements[:i] + elements[i+1:]):
            permutations.append([elements[i]] + p)
    return permutations

# Test functions and measure execution time
def test_function(func, *args):
    start_time = time.time()
    result = func(*args)
    end_time = time.time()
    print(f"Function: {func.__name__}\nResult: {result}\nTime taken: {end_time - start_time:.6f} seconds\n")

# Test data
elements = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5]
sorted_elements = sorted(elements)
target = 5

# Running the tests
test_function(is_first_element_null, elements)
test_function(find_max, elements)
test_function(bubble_sort, elements.copy())
test_function(binary_search, sorted_elements, target)
test_function(fibonacci, 10)  # Be cautious with larger numbers as it can take a very long time
test_function(permute, elements[:4])  # Limiting the size due to factorial time complexity