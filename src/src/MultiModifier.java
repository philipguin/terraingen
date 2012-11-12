package src;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MultiModifier<T> implements IModifier<T>
{
	private final List<IModifier<? super T>> modifiers = new ArrayList<IModifier<? super T>>();
	
	public MultiModifier(IModifier<? super T>... modifiers)
	{
		this.modifiers.addAll(Arrays.asList(modifiers));
	}
	
	public MultiModifier(List<? extends IModifier<? super T>> modifiers)
	{
		this.modifiers.addAll(modifiers);
	}
	
	public void modify(T values)
	{
		for (IModifier<? super T> modifier : modifiers)
			modifier.modify(values);
	}
}
