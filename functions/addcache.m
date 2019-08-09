function cache = addcache(cache,hashkey,input)

%Hard-coded limit to all function's caches
cachelimiter = 20;

%Check if current cache exceeds the number of cached states
if cache.size > cachelimiter
    
    %If so, get the hash table keys...
    keys = cache.keySet;
    keys = cell(keys.toArray);
    %... and remove the oldest cache entry 
    cache.remove(keys{end})
end

%Then add the new cache entry
cache.put(hashkey,input);


end