def main():
    import sys
    data = sys.stdin.read().splitlines()
    n = int(data[0])
    recipes = {}
    base_items = set()
    # Parse the recipes
    for i in range(1, 1+n):
        line = data[i].strip()
        if '=' in line:
            parts = line.split('=')
            result = parts[0]
            ingredients = parts[1].split('+')
            if result not in recipes:
                recipes[result] = []
            recipes[result].append(ingredients)
        # Also, all ingredients that are never results are base items.
    target = data[-1].strip()
    
    # We'll create a memo dictionary to store the minimum orbs for a potion.
    memo = {}
    
    def min_orbs(potion):
        if potion in memo:
            return memo[potion]
        # If the potion is not in recipes, it is a base item.
        if potion not in recipes:
            memo[potion] = 0
            return 0
        
        best = float('inf')
        for recipe in recipes[potion]:
            total_cost = 0
            valid = True
            for ing in recipe:
                cost_ing = min_orbs(ing)
                total_cost += cost_ing
            num_ing = len(recipe)
            total_cost += (num_ing - 1)
            if total_cost < best:
                best = total_cost
        memo[potion] = best
        return best
    
    result = min_orbs(target)
    print(result)

if __name__ == "__main__":
    main()