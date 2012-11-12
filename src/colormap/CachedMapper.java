package colormap;




public class CachedMapper<O> implements Mapper<Float, O>
{
	private final int range;
	private final O[] table;
	
	@SuppressWarnings("unchecked")
	public CachedMapper(Mapper<Float, ? extends O> toCache, int size)
	{
		this.range = size - 1;
		this.table = (O[])new Object[size];
		
		for (int i = 0; i < size; ++i)
			table[i] = toCache.map((float)i / range);
	}
	
	@Override
	public O map(Float value)
	{
		return table[(int)(value * range + .5f)];
	}
}
